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
getProjectSummary("TCGA-KIRC")

query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 100, directory="./TCGA/kirc")

data <- GDCprepare(query, directory = "./TCGA/kirc")
save(data, file = "data/kirc.Rdata")

load("data/kirc.Rdata")

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
  project = "TCGA-KIRC", 
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

save(all_data, file = "data/kirc_alldata.Rdata")


###
rm(list = ls())
load("data/kirc_alldata.Rdata")

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

###
rm(list = ls())
load("data/kirc_alldata.Rdata")

all_data$GADH11111 = (all_data$TSC22D1 + all_data$TMBIM6 + all_data$TM9SF2 + all_data$RBM47 + all_data$PPP1CB + 
                        all_data$HSD17B4 + all_data$`HLA-E` + all_data$ERGIC1 + all_data$CANX + all_data$ATP1B1) / 
                     (10 * all_data$ARGLU1)
cg = c("GADH11111")

data = all_data[-(1:4)]
data = log2(data + 1)
data = as.data.frame(t(data))

cg=cg[cg %in% rownames(data)]
cgl = list(cg=cg)

geneSet <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, cgl, names(cgl)))

data = as.matrix(data)
a = gsvaParam(
  data,
  geneSet,
  minSize = 1,
  maxSize = Inf,
  kcdf = c("Gaussian"),
  kcdfNoneMinSampleSize = 50,
  tau = 1,
  maxDiff = TRUE,
  absRanking = FALSE,
  sparse = TRUE
)

dat.gsva <- gsva(a)
dat.gsva=dat.gsva[1,]
group.gsva <- ifelse( dat.gsva < 0,'group1','group2')
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
life_data$futime = life_data$futime * 12

fit <- survfit(Surv(futime, fustat) ~ group, data = life_data)
pdf(file = "./Figure/step3_2.pdf", onefile = F, width = 8)
ggsurvplot(fit, data = life_data,
           size = 1,
           palette="lancet",
           conf.int = F,
           pval = TRUE,
           pval.method = TRUE,
           pval.size = 10,
           log.rank.weights = "1",
           legend.labs = c("NCG_low", "NCG_high"),
           xlim = c(0,120),
           break.time.by = 20,
           xlab = "Months",
           font.x = c(25, "bold", "black"),
           font.y = c(25, "bold", "black"),
           font.tickslab = c(20, "plain", "black"),
           legend.title = " ", 
           legend = "top",
           font.legend = 20,
           #title = sprintf("%s", gene),
           font.title = c(30, "bold", "black"),
           ggtheme = theme_classic2())
dev.off()

### immune analysis
rm(list = ls())

library(immunedeconv)
library(TCGAbiolinks)
library("tidyverse")
library(ggpubr)

load("data/kirc.Rdata")

tpm_ESEC = assay(data, "tpm_unstrand")

gene_id = rowData(data)
pd_data = colData(data)

gene_id_pc = gene_id %>% as_tibble() %>% filter(gene_type == "protein_coding")
tpm_ESEC = tpm_ESEC[gene_id_pc$gene_id,]

tpm_ESEC = as.data.frame(tpm_ESEC)

nor_data = tpm_ESEC  %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble() %>% 
  inner_join(gene_id_pc, by = "gene_id") %>% 
  as.data.frame()

nb = length(tpm_ESEC) + 1

nor_data$mean_expression <- rowMeans(nor_data[2:nb])
nor_data <- nor_data[order(nor_data$gene_name, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$gene_name), ]
rownames(nor_data) = nor_data$gene_name
nor_data = nor_data[2:nb]

delcol = which(is.na(nor_data), arr.ind = TRUE) %>% 
  as_tibble() %>% 
  dplyr::select(col) %>% 
  distinct(col)

#nor_data = nor_data[-delcol$col]


## load pd
query <- GDCquery(
  project = "TCGA-KIRC", 
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

save(all_data, file = "data/TPM_kirc_alldata.Rdata")

rm(list = ls())
load("data/TPM_kirc_alldata.Rdata")

data = all_data[-c(1:4)]
data = as.data.frame(t(data))

write.table(data, file = "data/KIRC_TPM.txt", sep = "\t")

LM22 = read.xlsx("./data/LM22.xls", "SuppTable1_GEP_Matrix")
rownames(LM22) = LM22$Gene.symbol
LM22 = LM22[-1]
write.table(LM22, file = "data/LM22.txt", sep = "\t")

source("../scRNA-snRNA/script/cibersort.R")

immune_cibersort <- CIBERSORT(
  sig_matrix = "data/LM22.txt", 
  mixture_file = "data/KIRC_TPM.txt",
  perm = 10, QN = F
)

save(immune_cibersort, file = "data/kirc_immune_cibersort.Rdata")

rm(list = ls())
load("data/kirc_immune_cibersort.Rdata")
load("data/kirc_alldata.Rdata")

all_data$GADH11111 = (all_data$TSC22D1 + all_data$TMBIM6 + all_data$TM9SF2 + all_data$RBM47 + all_data$PPP1CB + 
                        all_data$HSD17B4 + all_data$`HLA-E` + all_data$ERGIC1 + all_data$CANX + all_data$ATP1B1) / 
  (10 * all_data$ARGLU1)
cg = c("GADH11111")

data = all_data[-(1:4)]
data = log2(data + 1)
data = as.data.frame(t(data))

cg=cg[cg %in% rownames(data)]
cgl = list(cg=cg)

geneSet <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, cgl, names(cgl)))

data = as.matrix(data)
a = gsvaParam(
  data,
  geneSet,
  minSize = 1,
  maxSize = Inf,
  kcdf = c("Gaussian"),
  kcdfNoneMinSampleSize = 50,
  tau = 1,
  maxDiff = TRUE,
  absRanking = FALSE,
  sparse = TRUE
)

dat.gsva <- gsva(a)
dat.gsva=dat.gsva[1,]
group.gsva <- ifelse( dat.gsva < 0,"NCG_low", "NCG_high")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
life_data$group  = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

class(immune_cibersort)
immune_cibersort = as.data.frame(immune_cibersort)

new_data = immune_cibersort %>% 
  rownames_to_column("Barcode") %>% 
  as_tibble() %>% 
  inner_join(life_data, by = "Barcode") %>% 
  as.data.frame()

rownames(new_data) = new_data$Barcode
new_data = new_data[c(2:23, 30)]

y_positions <- new_data %>%
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>%
  group_by(Immune_Cells) %>%
  summarise(max_value = max(Proportion) * 1.3) %>% 
  pull(max_value)

pdf(file="./Figure/step3_10.pdf", width = 15)
new_data %>% 
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>% 
  ggplot(aes(x = Immune_Cells, y = Proportion)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format
dev.off()