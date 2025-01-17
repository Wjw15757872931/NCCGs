
## immune cycle  http://biocc.hrbmu.edu.cn/TIP/
rm(list = ls())

load("./data/tpm_step5_all_survial_data.Rdata")
data = all_data[-c(1:4)]
data = as.data.frame(t(data))

data_1 = data[1:50]
data_2 = data[51:100]
data_3 = data[101:150]
data_4 = data[151:200]
data_5 = data[201:233]

write.table(data_1, file = "data/CESC_TPM_1_50.txt", sep = "\t")
write.table(data_2, file = "data/CESC_TPM_51_100.txt", sep = "\t")
write.table(data_3, file = "data/CESC_TPM_101_150.txt", sep = "\t")
write.table(data_4, file = "data/CESC_TPM_151_200.txt", sep = "\t")
write.table(data_5, file = "data/CESC_TPM_201_233.txt", sep = "\t")

a1 = read.table("./data/ssGSEA.normalized.score1.txt", header = T, row.names = 1)
a2 = read.table("./data/ssGSEA.normalized.score2.txt", header = T, row.names = 1)
a3 = read.table("./data/ssGSEA.normalized.score3.txt", header = T, row.names = 1)
a4 = read.table("./data/ssGSEA.normalized.score4.txt", header = T, row.names = 1)
a5 = read.table("./data/ssGSEA.normalized.score5.txt", header = T, row.names = 1)

immune_cyc_score = cbind(cbind(cbind(a1, a2),cbind(a3, a4)), a5)

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

immune_cyc = as.data.frame(t(immune_cyc_score))

new_data = immune_cyc %>% 
  rownames_to_column("Barcode") %>% 
  as_tibble() %>% 
  mutate(Barcode = gsub("\\.", "-", Barcode)) %>% 
  inner_join(life_data, by = "Barcode") %>% 
  as.data.frame()

rownames(new_data) = new_data$Barcode
new_data = new_data[c(2:24, 28)]

y_positions <- new_data %>%
  pivot_longer(cols = -group, names_to = "immune_cyc_type", values_to = "Score") %>%
  group_by(immune_cyc_type) %>%
  summarise(max_value = max(Score) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("Step4.CD4_T_cell.recruiting", "Step4.CD8_T_cell.recruiting"
                  , "Step4.Monocyte.recruiting", "Step4.NK_cell.recruiting",
                  "Step4.Neutrophil.recruiting", "Step4.Th1_cell.recruiting",
                  "Step4.T_cell.recruiting", "Step5"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step3_13.pdf", width = 15)
tiff(file = "./Figure/step3_13.tiff", 
     width = 15 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
new_data %>% 
  pivot_longer(cols = -group, names_to = "immune_cyc_type", values_to = "Score") %>%
  ggplot(aes(x = immune_cyc_type, y = Score)) +
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
  ggtitle("TIP")
dev.off()

### Immune checkpoint and HLA gene
rm(list = ls())

immune_checkpoint <- c("ADORA2A", "BTLA", "BTNL2", "CD160", "CD200", "CD200R1", "CD244", "CD27", "CD274", "CD276", 
                       "CD28", "CD40", "CD40LG", "CD44", "CD48", "CD70", "CD80", "CD86", "CTLA4", "HAVCR2", "HHLA2", 
                       "ICOS", "ICOSLG", "IDO1", "IDO2", "KIR3DL1", "LAG3", "LAIR1", "LGALS9", "NELL1", "NRP1", 
                       "PDCD1", "PDCD1LG2", "SELPLG", "TIGIT", "TMIGD2", "TNFRSF14", "TNFRSF18", "TNFRSF25", 
                       "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF14", "TNFSF15", "TNFSF18", "TNFSF4", "TNFSF9", "VTCN1")

HLA_genes = c("HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", 
              "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", 
              "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB5", "HLA-DRB6", "HLA-E", 
              "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-L")

load("data/step5_all_survial_data.Rdata")
cg = c("TTC3", "MGST1", 'ACTR3')

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

data = all_data[-(1:4)]
data = log2(data + 1)
data = data[colnames(data) %in% immune_checkpoint]

data$group = as.character(group.gsva)
#data$group  = factor(data$group, levels = c("NCG_low", "NCG_high"))

y_positions <- data %>%
  pivot_longer(cols = -group, names_to = "genes", values_to = "expr") %>%
  group_by(genes) %>%
  summarise(max_value = max(expr) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("BTLA", "CD200", "CD276", "CD48", "CD70", "CD80", "CD86",
                  "CTLA4", "HAVCR2", "HHLA2", "ICOS", "IDO1", "IDO2", "LAG3", 
                  "LAIR1", "NRP1", "PDCD1", "TIGIT", "TNFRSF4"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step3_14.pdf", width = 20)
tiff(file = "./Figure/step3_14.tiff", 
     width = 15 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
data %>% 
  pivot_longer(cols = -group, names_to = "genes", values_to = "expr") %>%
  ggplot(aes(x = genes, y = expr)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format  +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("Immune Checkpoint")
dev.off()

### 
data = all_data[-(1:4)]
data = log2(data + 1)
data = data[colnames(data) %in% HLA_genes]

data$group = as.character(group.gsva)
#data$group  = factor(data$group, levels = c("NCG_low", "NCG_high"))

y_positions <- data %>%
  pivot_longer(cols = -group, names_to = "genes", values_to = "expr") %>%
  group_by(genes) %>%
  summarise(max_value = max(expr) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("HLA-B", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", 
                  "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", 
                  "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-F"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step3_15.pdf", width = 15)
tiff(file = "./Figure/step3_15.tiff", 
     width = 15 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
data %>% 
  pivot_longer(cols = -group, names_to = "genes", values_to = "expr") %>%
  ggplot(aes(x = genes, y = expr)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  plot.format  +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("Antigen Presentation")
dev.off()
