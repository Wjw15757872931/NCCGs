## 
rm(list = ls())

library(TCGAbiolinks)
library(tidyverse)
library(survival)
library(survminer)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(SummarizedExperiment)

getGDCprojects()$project_id 
getProjectSummary("TCGA-READ")

query <- GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 100, directory="./TCGA/read")

data <- GDCprepare(query, directory = "./TCGA/read")
save(data, file = "data/read.Rdata")

load("data/read.Rdata")

counts_ESEC = assay(data, "unstranded")

gene_id = rowData(data)
pd_data = colData(data)

gene_id_pc = gene_id %>% as_tibble() %>% filter(gene_type == "protein_coding")
counts_ESEC = counts_ESEC[gene_id_pc$gene_id,]

counts_ESEC = as.data.frame(counts_ESEC)

nor_data = counts_ESEC  %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble() %>% 
  inner_join(gene_id_pc, by = "gene_id") %>% 
  as.data.frame()

nb = length(counts_ESEC) + 1

nor_data$mean_expression <- rowMeans(nor_data[2:nb])
nor_data <- nor_data[order(nor_data$gene_name, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$gene_name), ]
rownames(nor_data) = nor_data$gene_name
nor_data = nor_data[2:nb]

delcol = which(is.na(nor_data), arr.ind = TRUE) %>% 
  as_tibble() %>% 
  dplyr::select(col) %>% 
  distinct(col)
#if delcol is NULL, not RUN
#nor_data = nor_data[-delcol$col]


## load pd
query <- GDCquery(
  project = "TCGA-READ", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR XML"
)

GDCdownload(query)
cli <- GDCprepare_clinic(query,'follow_up')
cli <- cli  %>% as_tibble() %>% 
  dplyr::select(bcr_followup_barcode,vital_status,                            
                days_to_death,days_to_last_followup) %>%
  distinct(bcr_followup_barcode, .keep_all = TRUE)

cli = cli %>% filter(vital_status %in% c("Alive", "Dead"))
dead_patient <-  cli %>%                                                        
  dplyr::filter(vital_status == 'Dead') %>%
  dplyr::select(-days_to_last_followup) %>%
  plyr::rename(c(bcr_followup_barcode = 'Barcode',
                 vital_status = 'fustat',
                 days_to_death='futime'
  )) %>%
  mutate(fustat=ifelse(fustat=='Dead',1,0))%>%
  mutate(futime=futime/365)
alive_patient <-  cli %>%                              
  dplyr::filter(vital_status == 'Alive') %>%
  dplyr::select(-days_to_death) %>%
  plyr::rename(c(bcr_followup_barcode = 'Barcode',
                 vital_status = 'fustat',
                 days_to_last_followup='futime'
  )) %>%
  mutate(fustat=ifelse(fustat=='Dead',1,0))%>%
  mutate(futime=futime/365)

survival_data <- rbind(dead_patient,alive_patient)

survival_data = na.omit(survival_data)
survival_data = as.data.frame(survival_data)
rownames(survival_data) = survival_data$Barcode

nor_data = t(nor_data)
nor_data = as.data.frame(nor_data)
nor_data$ID1 <- substr(rownames(nor_data), 1, 12)

survival_data$ID1 <- substr(rownames(survival_data), 1, 12)
survival_data  = survival_data %>% 
  arrange(ID1, desc(futime)) %>%
  distinct(ID1, .keep_all = T)
nor_data = nor_data %>% distinct(ID1, .keep_all = T)

all_data = survival_data %>%
  rownames_to_column() %>%
  as_tibble() %>%
  inner_join(nor_data, by = "ID1") %>%
  as.data.frame()
rownames(all_data) = all_data$rowname
all_data = all_data[-1]

save(all_data, file = "data/read_alldata.Rdata")


###
rm(list = ls())
load("data/read_alldata.Rdata")

data = all_data[-c(1:4)][colnames(all_data[-c(1:4)]) %in% sc_sn_genes$symbol]
data = as.data.frame(t(data))
data = log2(data + 1)

a = all_data$ID1
r <- length(row.names(data))
k <- length(data)
data$name <- row.names(data)
unicox <- data.frame(geneID = NA,HR = NA,downCI = NA,upCI = NA,P.val = NA)
s = 1
autogn = rownames(data)
sur = all_data[1:4] %>% rownames_to_column() %>% mutate(Sample = ID1)
sur$futime = sur$futime

# Single cox
while(r > 0){
  gene <- autogn[r]
  if (gene %in% row.names(data)){
    b <- data[data$name==gene,1:k]
    b <- as.numeric(b)
    ab <- data.frame(Sample = a,group = b)
    sur_r <- merge(sur,ab,by = "Sample")
    avg <- mean(b)
    sur_r$group[sur_r$group < avg] <- 0
    sur_r$group[sur_r$group >= avg] <- 1
    sur_r$group <- as.numeric(sur_r$group)
    res.cox <- coxph(Surv(futime, fustat) ~ group, data = sur_r)
    res <- summary(res.cox)
    unicox[s,1] <- gene
    unicox[s,2] <- res$conf.int[1]
    unicox[s,3] <- res$conf.int[3]
    unicox[s,4] <- res$conf.int[4]
    unicox[s,5] <- res$waldtest[3]
    s = s + 1
  }
  ;r <- r-1
}
unicox_sub <- unicox[unicox$P.val < 0.05,]
