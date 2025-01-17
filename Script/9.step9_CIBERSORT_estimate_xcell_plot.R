#### CIBERSORT

rm(list = ls())

library(immunedeconv)
library(TCGAbiolinks)
library("tidyverse")
library(ggpubr)
library(ggtext)

query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
)

GDCdownload(query)
cesc_data <- GDCprepare(query)
save(cesc_data, file = "cesc.Rdata")

load("./data/cesc.Rdata")
tpm_ESEC = assay(cesc_data, "tpm_unstrand")

gene_id = rowData(cesc_data)
pd_data = colData(cesc_data)

gene_id_pc = gene_id %>% as_tibble() %>% filter(gene_type == "protein_coding")
tpm_ESEC = tpm_ESEC[gene_id_pc$gene_id,]

tpm_ESEC = as.data.frame(tpm_ESEC)

nor_data = tpm_ESEC %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble() %>% 
  inner_join(gene_id_pc, by = "gene_id") %>% 
  as.data.frame()

colnames(nor_data)[308:316]

nor_data$mean_expression <- rowMeans(nor_data[2:310])
nor_data <- nor_data[order(nor_data$gene_name, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$gene_name), ]
rownames(nor_data) = nor_data$gene_name
nor_data = nor_data[2:310]

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
save(all_data, file="./data/tpm_step5_all_survial_data.Rdata")

load("./data/tpm_step5_all_survial_data.Rdata")
data = all_data[-c(1:4)]
data = as.data.frame(t(data))

write.table(data, file = "data/CESC_TPM.txt", sep = "\t")

LM22 = read.xlsx("./data/LM22.xls", "SuppTable1_GEP_Matrix")
rownames(LM22) = LM22$Gene.symbol
LM22 = LM22[-1]
write.table(LM22, file = "data/LM22.txt", sep = "\t")

source("../scRNA-snRNA/script/cibersort.R")

immune_cibersort <- CIBERSORT(
  sig_matrix = "data/LM22.txt", 
  mixture_file = "data/CESC_TPM.txt",
  perm = 10, QN = F
)

save(immune_cibersort, file = "data/immune_cibersort.Rdata")

load("data/immune_cibersort.Rdata")
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
  ifelse(x %in% c("NK.cells.resting", "T.cells.CD4.memory.activated", "T.cells.CD4.memory.resting", "T.cells.CD8"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=15),
                  axis.title=element_text(color="black",size=15),
                  plot.title=element_text(color="black",size=30, face = "bold"),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)


pdf(file="./Figure/step3_10.pdf", width = 15)
tiff(file = "./Figure/step3_10.tiff", 
     width = 15 * 300,  
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
esimate_score = read.table("./data/TCGA_CESC_ESTIMATE.txt", header = T)
load("./data/step5_all_survial_data.Rdata")

esimate_score$ID1 <- substr(esimate_score$ID, 1, 12)
pd_data = all_data[1:4]
n_data = esimate_score %>% 
  as_tibble() %>% 
  inner_join(all_data, by = "ID1") %>%
  distinct(Barcode, .keep_all = T) %>% 
  as.data.frame()

rownames(n_data) = n_data$Barcode

cg = c("TTC3", "MGST1", 'ACTR3')

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
                  plot.title=element_text(color="black",size=25, face = "bold"),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)

life_data = life_data[c(2:4, 9)]

life_data %>% 
  as_tibble() %>% 
  pivot_longer(cols = -group, names_to = "Type", values_to = "Score")

color_labels <- function(x) {
  ifelse(x %in% c("Immune_score", "ESTIMATE_score"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

y_positions <- life_data %>% 
  as_tibble() %>% 
  pivot_longer(cols = -group, names_to = "Type", values_to = "Score") %>%
  group_by(Type) %>%
  summarise(max_value = max(Score) * 1.3) %>% 
  pull(max_value)


pdf(file="./Figure/step3_11.pdf")
tiff(file = "./Figure/step3_11.tiff", 
     width = 6 * 300,   
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


a1 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_low") %>% 
  select(Immune_score)

a2 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_high") %>% 
  select(Immune_score)

wilcox.test(a1$Immune_score, a2$Immune_score)

pdf(file="./Figure/step3_7.pdf")
life_data[c(3,9)] %>% 
  ggplot(aes(x = group, y = Immune_score)) + 
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4)+
  ggtitle(label = expression(italic("P") == "0.0002465")) +
  plot.format +
  theme(legend.position = "none",
        axis.title.x = element_blank())
dev.off()

a1 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_low") %>% 
  select(Stromal_score)

a2 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_high") %>% 
  select(Stromal_score)

wilcox.test(a1$Stromal_score, a2$Stromal_score)

pdf(file="./Figure/step3_8.pdf")
life_data[c(2,9)] %>% 
  ggplot(aes(x = group, y = Stromal_score)) + 
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4)+
  ggtitle(label = expression(italic("P") == "0.9109")) +
  plot.format +
  theme(legend.position = "none",
        axis.title.x = element_blank())
dev.off()

a1 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_low") %>% 
  select(ESTIMATE_score)

a2 = life_data %>% 
  as_tibble() %>% 
  filter(group == "NCG_high") %>% 
  select(ESTIMATE_score)

wilcox.test(a1$ESTIMATE_score, a2$ESTIMATE_score)

pdf(file="./Figure/step3_9.pdf")
life_data[c(4,9)] %>% 
  ggplot(aes(x = group, y = ESTIMATE_score)) + 
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4)+
  ggtitle(label = expression(italic("P") == "0.01428")) +
  plot.format +
  theme(legend.position = "none",
        axis.title.x = element_blank())
dev.off()

## Xcell
rm(list = ls())
library(xCell)

load("./data/tpm_step5_all_survial_data.Rdata")
data = all_data[-c(1:4)]
data = as.data.frame(t(data))

xCell  = xCellAnalysis(data, 
                       rnaseq = TRUE, # RNA-seq
                       parallel.sz = 1, 
                       parallel.type = "SOCK") 


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
group.gsva <- ifelse( dat.gsva < 0,"LRG", "HRG")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
#life_data$group  = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

xCell = as.data.frame(t(xCell))

# 定义按名称相似性排序的免疫细胞顺序
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

y_positions <-  new_data%>%
  group_by(immune_cells) %>%
  summarise(max_value = max(immune_Score, na.rm = TRUE) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("B-cells", "Memory B-cells", "Class-switched memory B-cells", 
                  "CD4+ T-cells", "CD4+ naive T-cells", "CD4+ memory T-cells", "CD4+ Tem",
                  "CD8+ T-cells", "CD8+ naive T-cells", "CD8+ Tcm", "CD8+ Tem",
                  "Th1 cells",  "Tregs", "NKT",
                  "Macrophages", "Macrophages M1", "Macrophages M2",
                  "DC", "aDC", "cDC", "iDC", "pDC",
                  "Neutrophils", "Basophils",
                  "ImmuneScore", "MicroenvironmentScore"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step3_12.pdf", width = 20)
tiff(file = "./Figure/step3_12.tiff", 
     width = 20 * 300,  
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
  plot.format +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("Xcell")+
  labs(y = "Score")
dev.off()


