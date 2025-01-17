rm(list = ls())

library(TCGAbiolinks)
library(tidyverse)
library(survival)
library(survminer)
library(clusterProfiler)
library(GSVA)
library(GSEABase)

getGDCprojects()$project_id 
getProjectSummary("TCGA-HNSC")

query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 100, directory="./TCGA/hnsc")

data <- GDCprepare(query, directory = "./TCGA/hnsc")
save(data, file = "hnsc.Rdata")

load("data/hnsc.Rdata")

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

colnames(nor_data)[565:576]
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

nor_data = nor_data[-delcol$col]


## load pd
query <- GDCquery(
  project = "TCGA-HNSC", 
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

save(all_data, file = "data/hnsc_alldata.Rdata")

###
load("data/hnsc_alldata.Rdata")

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
tabletext <- cbind(
  c("GeneID", unicox_sub$geneID),
  c("HR (95% CI)", paste0(
    round(unicox_sub$HR, 2), " (",
    round(unicox_sub$downCI, 2), " - ",
    round(unicox_sub$upCI, 2), ")")
  ),
  c("P-value", formatC(unicox_sub$P.val, format = "e", digits = 2))
)

library(forestplot)

pdf(file = "./Figure/step7_2.pdf", onefile = F, width = 10)
forestplot(
  labeltext = tabletext,
  mean = c(NA, unicox_sub$HR),               
  lower = c(NA, unicox_sub$downCI),           
  upper = c(NA, unicox_sub$upCI),             
  zero = 1,                                 
  is.summary = c(TRUE, rep(FALSE, nrow(unicox_sub))), 
  xlog = TRUE,                        
  col = forestplot::fpColors(box = "black", line = "darkblue", summary = "royalblue"),
  xlab = "Hazard Ratio"
)
dev.off()


###
rm(list = ls())
load("data/hnsc_alldata.Rdata")

cg = c("TPT1", "TPM3", "TMBIM6", "TM9SF2", "PERP", "PDIA3", "PABPC1", 
       "P4HB", "GAPDH", "CANX", "CALR", "ARPC2", "ACTR3" )


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
group.gsva <- ifelse( dat.gsva < 0,'LRG','HRG')
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
life_data$group = factor(life_data$group, levels = c("LRG", "HRG"))
life_data$futime = life_data$futime * 12

cox_fit <- coxph(Surv(futime, fustat) ~ group, data = life_data)
hr <- round(summary(cox_fit)$conf.int[1],3)
hr_lower <- round(summary(cox_fit)$conf.int[3], 3)
hr_upper <- round(summary(cox_fit)$conf.int[4], 3)

pdf(file = "./Figure/step7_1.pdf", onefile = F, width = 8)
tiff(file = "./Figure/step7_1.tiff", 
     width = 10 * 300,   
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
survfit2(Surv(futime, fustat) ~ group, data = life_data) %>% 
  ggsurvfit(linewidth = 0.8) +
  add_pvalue(location  = "annotation", x = 180, y = 0.88 , hjust = 1,
             size = 4, caption = "Stratified log-rank {p.value}") +
  annotate("text",   
           x = 180,
           y = 0.95,
           label = paste("Hazard ratio,", round(hr,2), "(95% CI,", round(hr_lower,2), "-", round(hr_upper,2), ")", sep = ""),
           size = 4,
           hjust = 1
  ) +
  annotate("text",   
           x = 120,
           y = 0.5,
           label = "LRG",
           size = 4.5,
           hjust = 0,
  ) +  
  annotate("text",  
           x = 75,
           y = 0.25,
           label = "HRG",
           size = 4.5,
           hjust = 1  
  ) +  
  labs(title = "TCGA-HNSC(n=437)",
       x = "Months",
       y = "overall survival",
  ) +
  scale_x_continuous(breaks = seq(0, 200, 30), expand = c(0.03, 0)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  scale_color_manual(values = c('#2D5662', '#F7931D')) +
  scale_fill_manual(values = c('#2D5662', '#F7931D')) +
  theme_classic() +
  theme(axis.text = element_text(size = 14,color = "black"), 
        axis.title.y = element_text(size = 14,color = "black", vjust = 1),
        axis.title.x = element_text(size = 14,color = "black"),
        panel.grid.major.y =  element_line(color = "#DCDDDF"),
        legend.position = "none",
        plot.title=element_text(color="black",size=30, face = "bold")
        )
dev.off()

### immune analysis
rm(list = ls())

library(immunedeconv)
library(TCGAbiolinks)
library("tidyverse")
library(ggpubr)

load("data/hnsc.Rdata")

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

nor_data = nor_data[-delcol$col]

## load pd
query <- GDCquery(
  project = "TCGA-HNSC", 
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

save(all_data, file = "data/TPM_hnsc_alldata.Rdata")
load("data/TPM_hnsc_alldata.Rdata")

data = all_data[-c(1:4)]
data = as.data.frame(t(data))

write.table(data, file = "data/HNSC_TPM.txt", sep = "\t")

LM22 = read.xlsx("./data/LM22.xls", "SuppTable1_GEP_Matrix")
rownames(LM22) = LM22$Gene.symbol
LM22 = LM22[-1]
write.table(LM22, file = "data/LM22.txt", sep = "\t")

source("../scRNA-snRNA/script/cibersort.R")

immune_cibersort <- CIBERSORT(
  sig_matrix = "data/LM22.txt", 
  mixture_file = "data/HNSC_TPM.txt",
  perm = 10, QN = F
)

save(immune_cibersort, file = "data/hnsc_immune_cibersort.Rdata")

load("data/hnsc_immune_cibersort.Rdata")

cg = c("TPT1", "TPM3", "TMBIM6", "TM9SF2", "PERP", "PDIA3", "PABPC1", 
       "P4HB", "GAPDH", "CANX", "CALR", "ARPC2", "ACTR3" )

load("data/hnsc_alldata.Rdata")

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
group.gsva <- ifelse( dat.gsva < 0,"LRG", "HRG")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
#life_data$group  = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

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

color_labels <- function(x) {
  ifelse(x %in% c("B.cells.naive", "Dendritic.cells.activated", "Eosinophils",
                  "Mast.cells.resting", "Mast.cells.activated", "NK.cells.resting",
                  "Plasma.cells", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                  "T.cells.CD8", "T.cells.follicular.helper", "T.cells.gamma.delta",
                  "T.cells.regulatory..Tregs."), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step7_3.pdf", width = 15)
tiff(file = "./Figure/step7_3.tiff", 
     width = 20 * 300,   
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
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
  plot.format +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("CIBERSORT") +
  labs(y = "Score")
dev.off()

#### estimate_score
rm(list = ls())
esimate_score = read.table("./data/TCGA_HNSC_estimate.txt", header = T)
load("data/hnsc_alldata.Rdata")

esimate_score$ID1 <- substr(esimate_score$ID, 1, 12)
pd_data = all_data[1:4]
n_data = esimate_score %>% 
  as_tibble() %>% 
  inner_join(all_data, by = "ID1") %>%
  distinct(Barcode, .keep_all = T) %>% 
  as.data.frame()

rownames(n_data) = n_data$Barcode

cg = c("TPT1", "TPM3", "TMBIM6", "TM9SF2", "PERP", "PDIA3", "PABPC1", 
       "P4HB", "GAPDH", "CANX", "CALR", "ARPC2", "ACTR3" )

data = n_data[-(1:8)]
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
group.gsva <- ifelse(dat.gsva < 0, "LRG", "HRG")
table(group.gsva)

life_data = n_data[1:8]
life_data$group = as.character(group.gsva)
#life_data$group = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=15),
                  axis.title=element_text(color="black",size=15),
                  plot.title=element_text(color="black",size=15),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)

life_data = life_data[c(2:4, 9)]

life_data %>% 
  as_tibble() %>% 
  pivot_longer(cols = -group, names_to = "Type", values_to = "Score")

y_positions <- life_data %>% 
  as_tibble() %>% 
  pivot_longer(cols = -group, names_to = "Type", values_to = "Score") %>%
  group_by(Type) %>%
  summarise(max_value = max(Score) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("Immune_score", "ESTIMATE_score"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step7_4.pdf")
tiff(file = "./Figure/step7_4.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
life_data %>% 
  as_tibble() %>% 
  pivot_longer(cols = -group, names_to = "Type", values_to = "Score") %>%
  ggplot(aes(x = Type, y = Score)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format+
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("ESTIMATE")
dev.off()

## Xcell
rm(list = ls())

library(xCell)

load("data/TPM_hnsc_alldata.Rdata")
data = all_data[-c(1:4)]
data = as.data.frame(t(data))

xCell  = xCellAnalysis(data, 
                       rnaseq = TRUE, 
                       parallel.sz = 1, 
                       parallel.type = "SOCK") 

cg = c("TPT1", "TPM3", "TMBIM6", "TM9SF2", "PERP", "PDIA3", "PABPC1", 
       "P4HB", "GAPDH", "CANX", "CALR", "ARPC2", "ACTR3" )

load("data/hnsc_alldata.Rdata")

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
group.gsva <- ifelse( dat.gsva < 0,"LRG", "HRG")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
#life_data$group  = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

xCell = as.data.frame(t(xCell))

ordered_immune_cells <- c("B-cells", "Memory B-cells", "naive B-cells", "Class-switched memory B-cells", 
                          "CD4+ T-cells", "CD4+ naive T-cells", "CD4+ memory T-cells", "CD4+ Tcm", "CD4+ Tem",
                          "CD8+ T-cells", "CD8+ naive T-cells", "CD8+ Tcm", "CD8+ Tem",
                          "Th1 cells", "Th2 cells", "Tregs", "Tgd cells",
                          "NK cells", "NKT",
                          "Macrophages", "Macrophages M1", "Macrophages M2",
                          "DC", "aDC", "cDC", "iDC", "pDC",
                          "Monocytes", "Neutrophils", "Basophils", "Eosinophils", "Mast cells",
                          "ImmuneScore", "StromaScore", "MicroenvironmentScore", "group")

new_data = xCell %>% 
  rownames_to_column("Barcode") %>% 
  as_tibble() %>% 
  inner_join(life_data, by = "Barcode") %>% 
  dplyr::select(all_of(ordered_immune_cells)) %>% 
  pivot_longer(-group, names_to = "immune_cells", values_to = "immune_Score")


new_data$immune_cells = factor(new_data$immune_cells, levels = ordered_immune_cells[-36])

color_labels <- function(x) {
  ifelse(x %in% c("B-cells", "Memory B-cells", "naive B-cells", "Class-switched memory B-cells", 
                  "CD4+ T-cells", "CD4+ naive T-cells", "CD4+ memory T-cells",
                  "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem",
                  "Th1 cells", "Th2 cells", "Tgd cells",
                  "aDC", "pDC",
                  "Monocytes", "Basophils", "Mast cells",
                  "ImmuneScore", "StromaScore", "MicroenvironmentScore"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

y_positions <-  new_data%>%
  group_by(immune_cells) %>%
  summarise(max_value = max(immune_Score, na.rm = TRUE) * 1.3) %>% 
  pull(max_value)

pdf(file="./Figure/step7_5.pdf", width = 20)
tiff(file = "./Figure/step7_5.tiff", 
     width = 25 * 300,   
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
new_data %>% 
  ggplot(aes(x = immune_cells, y = immune_Score)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format+
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("Xcell")+
  labs(y = "Score")
dev.off()
