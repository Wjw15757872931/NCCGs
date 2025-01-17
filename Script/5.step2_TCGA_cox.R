## TCGA
rm(list = ls())

library(TCGAbiolinks)
library("tidyverse")
library(SummarizedExperiment)
library("survival")
library("survminer")
library("clusterProfiler")

getGDCprojects()$project_id 
getProjectSummary("TCGA-CESC")

query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 100, directory="./bioproject/TCGA/RNA-seq/raw")

data <- GDCprepare(query, save = TRUE, directory = "./bioproject/TCGA/RNA-seq/raw", save.filename = "CESC_hg38.RData")
counts_ESEC = assay(data, "unstranded")

gene_id = rowData(data)
pd_data = colData(data)

gene_id_pc = gene_id %>% as_tibble() %>% filter(gene_type == "protein_coding")
counts_ESEC = counts_ESEC[gene_id_pc$gene_id,]

which(pd_data$tissue_type == "Normal")
which(pd_data$tissue_type == "Tumor")

Normal_count_data = counts_ESEC[,rownames(pd_data[which(pd_data$tissue_type == "Normal"),])]
Tumor_count_data = counts_ESEC[,rownames(pd_data[which(pd_data$tissue_type == "Tumor"),])]
save(Normal_count_data, Tumor_count_data, gene_id_pc, file = "step2.Rdata")

load("data/step2.Rdata")
rm(Normal_count_data)

Tumor_count_data[1:6, 1:6]
gene_id_pc[1:6,1:6]
Tumor_count_data = as.data.frame(Tumor_count_data)
#nor_data = log2(Tumor_count_data + 1)
nor_data = Tumor_count_data %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble() %>% 
  inner_join(gene_id_pc, by = "gene_id") %>% 
  as.data.frame()

nb = length(Tumor_count_data) + 1

nor_data$mean_expression <- rowMeans(nor_data[2:nb])
nor_data <- nor_data[order(nor_data$gene_name, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$gene_name), ]
rownames(nor_data) = nor_data$gene_name
nor_data = nor_data[2:nb]


## load pd
query <- GDCquery(
  project = "TCGA-CESC", 
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
save(survival_data, file = "TCGA_survival_data.Rdata")

################################
load("data/TCGA_survival_data.Rdata")
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
save(all_data, file="./data/step5_all_survial_data.Rdata")

load("./data/step5_all_survial_data.Rdata")
data = all_data[-(1:4)][colnames(all_data[-(1:4)]) %in% sc_sn_genes$symbol]
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
# 绘制森林图

pdf(file = "./Figure/step3_1.pdf", onefile = F, width = 10)
tiff(file = "./Figure/step3_1.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
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

####################################################
## live analysis
load("./data/step5_all_survial_data.Rdata")

library("survival")
library("survminer")
library("clusterProfiler")
library(tidyverse)

head(all_data)
hub_gene = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1")

data = all_data[-(1:4)]
data = log2(data + 1)
data = as.data.frame(t(data))

a = all_data$ID1
r <- length(hub_gene)
k <- length(data)
data$name <- row.names(data)
sur = all_data[1:4] %>% rownames_to_column() %>% mutate(Sample = ID1)
sur$futime = sur$futime * 12

while (r > 0) {
  gene <- hub_gene[r]
  if (gene %in% row.names(data)){
    b <- data[data$name==gene,1:k]
    b <- as.numeric(b)
    ab <- data.frame(Sample = a,group = b)
    sur_r <- merge(sur,ab,by = "Sample")
    avg <- mean(b)
    sur_r$group[sur_r$group < avg] <- 0
    sur_r$group[sur_r$group >= avg] <- 1
    sur_r$group <- as.numeric(sur_r$group)
    fit <- survfit(Surv(futime, fustat) ~ group, data = sur_r)
    pdf(sprintf("live plot %s.pdf", gene), width = 8, onefile = FALSE)
    p = ggsurvplot(fit, data = sur_r,
                   size = 1,
                   palette="lancet",
                   conf.int = F,
                   pval = TRUE,
                   pval.method = TRUE,
                   pval.size = 10,
                   log.rank.weights = "1",
                   legend.labs = c("low", "high"),
                   xlim = c(0,120),
                   break.time.by = 20,
                   xlab = "Months",
                   font.x = c(25, "bold", "black"),
                   font.y = c(25, "bold", "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend.title = " ", 
                   legend = "top",
                   font.legend = 20,
                   title = sprintf("%s", gene),
                   font.title = c(30, "bold", "black"),
                   ggtheme = theme_classic2())
    print(p)
    dev.off()
    r = r - 1
  }
}

### GSVA
rm(list = ls())

library("survival")
library("survminer")
library("clusterProfiler")
library(tidyverse)
library(GSVA)
library(GSEABase)

cg = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1")
cg = c("TTC3", "MGST1", 'ACTR3')

load("data/step5_all_survial_data.Rdata")

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
life_data$group = factor(life_data$group, levels = c('LRG','HRG'))
life_data$futime = life_data$futime * 12

cox_fit <- coxph(Surv(futime, fustat) ~ group, data = life_data)
hr <- round(summary(cox_fit)$conf.int[1],3)
hr_lower <- round(summary(cox_fit)$conf.int[3], 3)
hr_upper <- round(summary(cox_fit)$conf.int[4], 3)

pdf(file = "./Figure/step7_10.pdf", onefile = F, width = 8)
pdf(file = "./Figure/step3_2.pdf", onefile = F, width = 8)
tiff(file = "./Figure/step3_2.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
tiff(file = "./Figure/step7_10.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
survfit2(Surv(futime, fustat) ~ group, data = life_data) %>% 
  ggsurvfit(linewidth = 0.8) +
  add_pvalue(location  = "annotation", x = 200, y = 0.88 , hjust = 1,
             size = 4, caption = "Stratified log-rank {p.value}") +
  annotate("text",    
           x = 200,
           y = 0.95,
           label = paste("Hazard ratio,", round(hr,2), "(95% CI,", round(hr_lower,2), "-", round(hr_upper,2), ")", sep = ""),
           size = 4,
           hjust = 1
  ) +
  annotate("text",    
           x = 120,
           y = 0.78,
           label = "LRG",
           size = 4.5,
           hjust = 0,
  ) +  
  annotate("text",  
           x = 60,
           y = 0.55,
           label = "HRG",
           size = 4.5,
           hjust = 1  
  ) +  
  labs(title = "TCGA-CESC(n=233)",
       x = "Months",
       y = "overall survival",
  ) +
  scale_x_continuous(breaks = seq(0, 220, 30), expand = c(0.03, 0)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  scale_color_manual(values = c('#2D5662', '#F7931D')) +
  scale_fill_manual(values = c('#2D5662', '#F7931D')) +
  theme_classic() +
  theme(axis.text = element_text(size = 14,color = "black"), 
        axis.title.y = element_text(size = 14,color = "black", vjust = 1),
        axis.title.x = element_text(size = 14,color = "black"),
        panel.grid.major.y =  element_line(color = "#DCDDDF"),
        legend.position = "none",
        plot.title = element_text(size = 30, face = "bold"
        ))
dev.off()

###
fit <- survfit(Surv(futime, fustat) ~ group, data = life_data)
#pdf(file = "./Figure/step3_2.pdf", onefile = F, width = 8)
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
###
rm(list = ls())
library(xlsx)
library(tidyverse)
library(ggtext)

load("data/step6_all_data.Rdata")

data = all_data[-(1:5)][colnames(all_data[-(1:5)]) %in% sc_sn_genes$symbol]
data = log2(data + 1)
data$group = all_data$group
data$group = ifelse(data$group == "NCG_low", "LRG", "HRG")

y_positions <- data %>%
  pivot_longer(cols = -group, names_to = "genes", values_to = "Value") %>%
  group_by(genes) %>%
  summarise(max_value = max(Value) * 1.3) %>% 
  pull(max_value)


color_labels <- function(x) {
  ifelse(x %in% c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step6_2.pdf", width = 20)

data %>% 
  pivot_longer(cols = -group, names_to = "genes", values_to = "Value") %>%
  ggplot(aes(x = genes, y = Value)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y = y_positions) + # 使用更高的 y 位置
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("TCGA-CESC(n=233)")

dev.off()

#### HPV+ TCGA
rm(list = ls())
library(xlsx)
library(tidyverse)

load("data/step5_all_survial_data.Rdata")

hpvd = read.xlsx("F:/bioinformatics/f5.xlsx", "HPV main")
hpvd = na.omit(hpvd[1:7])
hpvd = hpvd[-which(hpvd$HPV.types == "negative"),]

all_data = all_data %>% 
  as_tibble() %>% 
  filter(ID1 %in% hpvd$Sample_ID)

### GSVA
library("survival")
library("survminer")
library("clusterProfiler")
library(tidyverse)
library(GSVA)
library(GSEABase)
library(ggsurvfit)

cg = c("TTC3", "ACTR3", "MGST1")
all_data$GADH11111 = ((all_data$TTC3 + all_data$DNAJC3 + all_data$ACTR3 + all_data$MGST1)) / (4*(all_data$TMSB4X))
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
group.gsva <- ifelse( dat.gsva < 0,'LRG','HRG')
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
life_data$group = factor(life_data$group, levels = c('LRG','HRG'))
life_data$futime = life_data$futime * 12

cox_fit <- coxph(Surv(futime, fustat) ~ group, data = life_data)
hr <- round(summary(cox_fit)$conf.int[1],3)
hr_lower <- round(summary(cox_fit)$conf.int[3], 3)
hr_upper <- round(summary(cox_fit)$conf.int[4], 3)

pdf(file = "./Figure/step7_11_hpv_all.pdf", onefile = F, width = 8)
pdf(file = "./Figure/step3_4_hpv.pdf", onefile = F, width = 8)
tiff(file = "./Figure/step3_4_hpv.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
tiff(file = "./Figure/step7_11_hpv.tiff", 
     width = 10 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
survfit2(Surv(futime, fustat) ~ group, data = life_data) %>% 
  ggsurvfit(linewidth = 0.8) +
  add_pvalue(location  = "annotation", x = 200, y = 0.88 , hjust = 1,
             size = 4, caption = "Stratified log-rank {p.value}") +
  annotate("text",    
           x = 200,
           y = 0.95,
           label = paste("Hazard ratio,", round(hr,2), "(95% CI,", round(hr_lower,2), "-", round(hr_upper,2), ")", sep = ""),
           size = 4,
           hjust = 1
           ) +
  annotate("text",    
           x = 120,
           y = 0.75,
           label = "LRG",
           size = 4.5,
           hjust = 0,
           ) +  
  annotate("text",   
           x = 70,
           y = 0.5,
           label = "HRG",
           size = 4.5,
           hjust = 1  
           ) +  
  labs(title = "HPV+ TCGA-CESC(n=153)",
       x = "Months",
       y = "overall survival",
       ) +
  scale_x_continuous(breaks = seq(0, 220, 30), expand = c(0.03, 0)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  scale_color_manual(values = c('#2D5662', '#F7931D')) +
  scale_fill_manual(values = c('#2D5662', '#F7931D')) +
  theme_classic() +
  theme(axis.text = element_text(size = 14,color = "black"), 
        axis.title.y = element_text(size = 14,color = "black", vjust = 1),
        axis.title.x = element_text(size = 14,color = "black"),
        panel.grid.major.y =  element_line(color = "#DCDDDF"),
        legend.position = "none",
        plot.title = element_text(size = 30, face = "bold") 
        )
dev.off()


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
# 绘制森林图

pdf(file = "./Figure/step3_5_hpv.pdf", onefile = F, width = 10)
tiff(file = "./Figure/step3_5_hpv.tiff", 
     width = 10 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
forestplot(
  labeltext = tabletext,
  mean = c(NA, unicox_sub$HR),              
  lower = c(NA, unicox_sub$downCI),           
  upper = c(NA, unicox_sub$upCI),            
  zero = 1,                                  
  is.summary = c(TRUE, rep(FALSE, nrow(unicox_sub))), 
  xlog = TRUE,                                
  col = forestplot::fpColors(box = "black", line = "darkblue", summary = "royalblue"),
  xlab = "Hazard Ratio",
  #title = "HPV+ TCGA-CESC(n=153)",                  # 添加标题
  #txt_gp = forestplot::fpTxtGp(
  #  title = grid::gpar(fontsize = 30, fontface = "bold", just = "right")
  #)  
)
dev.off()

## live analysis
load("./data/step5_all_survial_data.Rdata")

library("survival")
library("survminer")
library("clusterProfiler")
library(tidyverse)

head(all_data)
hub_gene = c("TTC3", "TMSB4X", "DNAJC3", "ACTR3", "MGST1")

data = all_data[-(1:4)]
data = log2(data + 1)
data = as.data.frame(t(data))

a = all_data$ID1
r <- length(hub_gene)
k <- length(data)
data$name <- row.names(data)
sur = all_data[1:4] %>% rownames_to_column() %>% mutate(Sample = ID1)
sur$futime = sur$futime * 12

while (r > 0) {
  gene <- hub_gene[r]
  if (gene %in% row.names(data)){
    b <- data[data$name==gene,1:k]
    b <- as.numeric(b)
    ab <- data.frame(Sample = a,group = b)
    sur_r <- merge(sur,ab,by = "Sample")
    avg <- mean(b)
    sur_r$group[sur_r$group < avg] <- 0
    sur_r$group[sur_r$group >= avg] <- 1
    sur_r$group <- as.numeric(sur_r$group)
    fit <- survfit(Surv(futime, fustat) ~ group, data = sur_r)
    pdf(sprintf("live plot hpv %s.pdf", gene), width = 8, onefile = FALSE)
    p = ggsurvplot(fit, data = sur_r,
                   size = 1,
                   palette="lancet",
                   conf.int = F,
                   pval = TRUE,
                   pval.method = TRUE,
                   pval.size = 10,
                   log.rank.weights = "1",
                   legend.labs = c("low", "high"),
                   xlim = c(0,120),
                   break.time.by = 20,
                   xlab = "Months",
                   font.x = c(25, "bold", "black"),
                   font.y = c(25, "bold", "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend.title = " ", 
                   legend = "top",
                   font.legend = 20,
                   title = sprintf("%s", gene),
                   font.title = c(30, "bold", "black"),
                   ggtheme = theme_classic2())
    print(p)
    dev.off()
    r = r - 1
  }
}

#####
rm(list = ls())
library(xlsx)
library(tidyverse)
library(ggtext)

load("data/step6_all_data.Rdata")
hpvd = read.xlsx("F:/bioinformatics/f5.xlsx", "HPV main")
hpvd = na.omit(hpvd[1:7])
hpvd = hpvd[-which(hpvd$HPV.types == "negative"),]

all_data = all_data %>% 
  as_tibble() %>% 
  filter(ID1 %in% hpvd$Sample_ID)

data = all_data[-(1:5)][colnames(all_data[-(1:5)]) %in% sc_sn_genes$symbol]
data = log2(data + 1)
data$group = all_data$group
data$group = ifelse(data$group == "NCG_low", "LRG", "HRG")

y_positions <- data %>%
  pivot_longer(cols = -group, names_to = "genes", values_to = "Value") %>%
  group_by(genes) %>%
  summarise(max_value = max(Value) * 1.3) %>% 
  pull(max_value)


color_labels <- function(x) {
  ifelse(x %in% c('TTC3', 'TMSB4X', 'MGST1', 'DNAJC3', 'ACTR3', 'ZFP36L1', 'TPT1', 'TPM3', 'TM9SF2', 'P4HB'), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step6_1_hpv.pdf", width = 20)
tiff(file = "./Figure/step6_1_hl.tiff", 
     width = 20 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")

data %>% 
  pivot_longer(cols = -group, names_to = "genes", values_to = "Value") %>%
  ggplot(aes(x = genes, y = Value)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y = y_positions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("TCGA-CESC(n=233)")

dev.off()
