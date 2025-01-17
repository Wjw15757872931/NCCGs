## inferCNV

rm(list = ls())

library(infercnv)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(BiocParallel)
library(AnnoProbe)

options(scipen = 100)
register(MulticoreParam(workers = 8, progressbar = TRUE))
memory.limit(size = size_in_megabytes * 1024^2)

sce_harmony = readRDS("step1data.rds")

epithelial_cells <- subset(sce_harmony, CellType == "Epithelial cells")

pc.num = 1:12
epithelial_cells <- RunUMAP(epithelial_cells, reduction = "harmony", dims = pc.num)
epithelial_cells <- RunTSNE(epithelial_cells, reduction = "harmony", dims = pc.num)
epithelial_cells <- FindNeighbors(epithelial_cells, reduction = "harmony", dims = pc.num)
epithelial_cells <- FindClusters(epithelial_cells, reduction = "harmony", resolution = 0.2)
saveRDS(epithelial_cells, "step2data.rds")

new_epithelial_clusters <- Idents(epithelial_cells)
sce_harmony$CellType_refined <- sce_harmony$CellType
sce_harmony$CellType_refined[names(new_epithelial_clusters)] <- as.character(new_epithelial_clusters)

sce_harmony = subset(sce_harmony, CellType %in% c("lymphocytes","Epithelial cells"))
sce_harmony$CellType_refined = factor(sce_harmony$CellType_refined, levels = c("0", "1", "2",
                                                                               "3", "4", "5",
                                                                               "6", "7", "lymphocytes",
                                                                               "Endothelial cells"))

dfcount = GetAssayData(sce_harmony, slot="data")
#dfcount = dfcount[VariableFeatures(sce_harmony),]
#dfcount = round(dfcount, digits=3)

groupinfo= data.frame(cellId = colnames(dfcount),cellType= sce_harmony$CellType_refined )

geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]

colnames(geneInfor) = c("V1", "V2", "V3", "V4")
geneInfor$V5 = str_extract(geneInfor$V2, "\\d+")
geneInfor$V5 = as.numeric(geneInfor$V5)
geneInfor = arrange(geneInfor, V5, V3)
geneInfor = geneInfor %>% as_tibble() %>% filter(V2 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")) %>% as.data.frame()
geneInfor = geneInfor[-5]

dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match(geneInfor[,1], rownames(dfcount) ),] 

#write.table(dfcount ,file ='data/expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = 'data/metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file ='data/geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)

rm(sce_harmony, groupinfo, geneInfor)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dfcount,
                                    annotations_file="data/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= 'data/geneFile.txt',
                                    ref_group_names="lymphocytes")

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, #10x 0.1 
                             out_dir='data', 
                             cluster_by_groups=TRUE,
                             analysis_mode="subclusters",
                             denoise=T,
                             HMM=T,
                             output_format = "pdf",
                             num_threads=8)

##########################################################
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

infercnv_obj = readRDS("./data/run.final.infercnv_obj")

expr <- infercnv_obj@expr.data

Group <- character(length(colnames(expr)))
Group[infercnv_obj@reference_grouped_cell_indices$lymphocytes] <- "lymphocytes"
Group[infercnv_obj@observation_grouped_cell_indices$`0`] <- "0"
Group[infercnv_obj@observation_grouped_cell_indices$`1`] <- "1"
Group[infercnv_obj@observation_grouped_cell_indices$`2`] <- "2"
Group[infercnv_obj@observation_grouped_cell_indices$`3`] <- "3"
Group[infercnv_obj@observation_grouped_cell_indices$`4`] <- "4"
Group[infercnv_obj@observation_grouped_cell_indices$`5`] <- "5"
Group[infercnv_obj@observation_grouped_cell_indices$`6`] <- "6"
Group[infercnv_obj@observation_grouped_cell_indices$`7`] <- "7"


rm(infercnv_obj)
CNV_data <- data.frame(cells = colnames(expr), Group = Group)

CNV_data$CNV_score= colMeans((expr)^2)
rm(infercnv_obj, expr, Group)

color_v = c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
            '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
            '#d69971','#df5734')

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

load("./data/CNV_DATA.Rdata")
pdf(file="Figure/step2_1.pdf", width = 10)
CNV_data %>% ggplot(aes(Group,CNV_score))+geom_violin(aes(fill=Group),color="NA")+
  scale_fill_manual(values = color_v) +
  plot.format +
  ylim(0.9, 1.2)
dev.off()

####
##sccancer
library(Seurat)

sce_harmony = readRDS("step1data.rds")
epithelial_cells = readRDS("step2data.rds")
source("sccancer.R")

new_epithelial_clusters <- Idents(epithelial_cells)
sce_harmony$CellType_refined <- sce_harmony$CellType
sce_harmony$CellType_refined[names(new_epithelial_clusters)] <- as.character(new_epithelial_clusters)

sce_harmony = subset(sce_harmony, CellType %in% c("lymphocytes","Epithelial cells"))
sce_harmony$CellType_refined = factor(sce_harmony$CellType_refined, levels = c("0", "1", "2",
                                                                               "3", "4", "5",
                                                                               "6", "7", "lymphocytes"
                                                                               ))

sccancer_score = runStemness(sce_harmony@assays$RNA@data)

save(sccancer_score, file="sccancer_stemness.Rdata")
load("./data/sccancer_stemness.Rdata")
CNV_data$sccancer_score = sccancer_score

pdf(file="Figure/step2_2.pdf", width = 10)
CNV_data %>% ggplot(aes(Group,sccancer_score))+geom_violin(aes(fill=Group),color="NA")+
  scale_fill_manual(values = color_v) +
  plot.format
dev.off()

####
library(ggsci)
library(Seurat)
library(tidyverse)
library(SCP)
library(BiocParallel)

register(MulticoreParam(workers = 8, progressbar = TRUE))
memory.limit(size = size_in_megabytes * 1024^2)

col <- pal_aaas()(9)
zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')

sce_harmony = readRDS("step1data.rds")
epithelial_cells = readRDS("step2data.rds")

new_epithelial_clusters <- Idents(epithelial_cells)
sce_harmony$CellType[names(new_epithelial_clusters)] <- as.character(new_epithelial_clusters)
sce_harmony$CellType = factor(sce_harmony$CellType, levels = c("lymphocytes", "0", "1",
                                                               "2", "3", "4", "5",
                                                               "6", "7", "Myeloid cells",
                                                               "Fibroblasts", "Endothelial cells", "Smooth muscle cells"
))

pdf(file="./Figure/step2_3.pdf")
CellDimPlot(                                                                           
  srt = sce_harmony, group.by = c("CellType"),
  reduction = "tsne", palcolor = zzm60colors                    
)
dev.off()


epithelial_cells$Cluster = epithelial_cells$seurat_clusters
epithelial_cells$Cluster = factor(epithelial_cells$Cluster, levels = c("0", "1",
                                                                       "2", "3",
                                                                       "4", "5",
                                                                       "6", "7"))
pdf(file="./Figure/step2_4.pdf")
CellDimPlot(                                                                           
  srt = epithelial_cells, group.by = c("CellType"),
  reduction = "tsne", palcolor = zzm60colors[2:9]                    
)
dev.off()

####
epi$Cluster = epi$seurat_clusters
Cellratio = prop.table(table(epi$Cluster, epi$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
saveRDS(Cellratio, file="./data/epi_sample_cellpro.rds")

Cellratio = readRDS("./Rdata/epi_sample_cellpro.rds")
Cellratio = as.data.frame(Cellratio)
colnames(Cellratio) = c("Cluster", "Group", "Freq")

pdf(file = "./Figure/step2_5.pdf", width = 10, height = 8)
ggplot(Cellratio,aes(Group, Freq,fill = Cluster)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cluster",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  #coord_flip() +
  theme(legend.position = 'right') +
  plot.format + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = zzm60colors[2:9])
dev.off()

####
epithelial_cells = readRDS("step2data.rds")
epithelial_cells$Cluster = epithelial_cells$seurat_clusters
epithelial_cells$Cluster = factor(epithelial_cells$Cluster, levels = c("0", "1",
                                                                       "2", "3",
                                                                       "4", "5",
                                                                       "6", "7"))
epi_degs <- FindAllMarkers(epithelial_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(epi_degs, file = "epi_deg.Rdata")

pdf(file = "step2_6.pdf", width = 10, onefile = F)
FeatureStatPlot(epithelial_cells,
                stat.by = c(
                  "KRT6A", 'KRT15', 'KRT17', 'CDKN2A', 'ERBB2',
                  "KRT5", "KRT8", 'TP63',
                  'MKI67', 'TOP2A'
                ),
                legend.position = "right", legend.direction = "vertical",
                group.by = "Cluster", stack = TRUE,
                palcolor = zzm60colors[2:9], ylab = "", xlab = ""
                
)
dev.off()

####
library(Seurat)
library(tidyverse)
library(SCP)

epi = readRDS("./data/step2data.rds")
table(epi$seurat_clusters, epi$Group)

VlnPlot(epi,
        features = c("KRT6A", 'KRT15', 'KRT17', 'CDKN2A', 'ERBB2'), 
        #features = c("KRT5", "KRT8", 'TP63'), 
        #features = c('MKI67', 'TOP2A'),
        split.by = 'seurat_clusters',stack = TRUE, 
        flip = T) + 
  theme(legend.position='none')+labs(x='',y='')

genes = c("KRT6A", 'KRT15', 'KRT17', 'CDKN2A', 'ERBB2',
          "KRT5", "KRT8", 'TP63',
          'MKI67', 'TOP2A')

load("./data/epi_deg.Rdata")
top3 = epi_degs %>% as_tibble() %>% 
  group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

plots <- lapply(top3$gene, function(gene) {
  FeaturePlot(epi, features = gene, reduction = "tsne" ,pt.size = 1) + 
    ggtitle(gene) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.size = unit(0.2, "cm"), 
          legend.text = element_text(size = 4),
          axis.title=element_text(color="black",size=7),
          plot.title=element_text(color="black",size=7),
    )
  
})

combined_plot <- wrap_plots(plots, ncol = 6)

pdf("./Figure/step2_7.pdf", width = 13)
print(combined_plot)
dev.off()

### gsva HALL-MARK GENES

library(Seurat)
library(GSVA)
library(irGSEA)
library(SeuratData)
library(RcppML)
library(irGSEA)

epi = readRDS("./data/step2data.rds")
load("./data/canceSEAgene_sets.Rdata")

epi <- irGSEA.score(object = epi,assay = "RNA",
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T,
                             species = "Homo sapiens", category = "H",
                             subcategor = NULL, geneid = "symbol",
                             method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                             kcdf = 'Gaussian')

epi$Cluster = epi$seurat_clusters
result.dge <- irGSEA.integrate(object = epi,
                               group.by = "Cluster",
                               method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"))

load("data/result.dge.Rdata")
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>%
  dplyr::pull(Name) %>% unique(.)

pdf("./Figure/step2_10.pdf", width = 10)
irGSEA.heatmap(object = result.dge,
                        method = "RRA",
                        top = 20,
                        show.geneset = NULL,
                        cluster.color = zzm60colors[2:9])
dev.off()

###########################################################
library(SCP)
library(tidyverse)
library(Seurat)

zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')

sce_harmony = readRDS("step1data.rds")
epithelial_cells = readRDS("step2data.rds")

new_epithelial_clusters <- Idents(epithelial_cells)
sce_harmony$CellType_refined <- sce_harmony$CellType
sce_harmony$CellType_refined[names(new_epithelial_clusters)] <- as.character(new_epithelial_clusters)

sce_harmony = subset(sce_harmony, CellType %in% c("lymphocytes","Epithelial cells"))
sce_harmony$CellType_refined = factor(sce_harmony$CellType_refined, levels = c("lymphocytes", "0", "1", "2",
                                                                               "3", "4", "5", "6", "7", 
))
sce_harmony$CellType = sce_harmony$CellType_refined

load("./data/CNV_DATA.Rdata")
sce_harmony$CNV_score = CNV_data$CNV_score

pdf("step2_8.pdf")
FeatureStatPlot(sce_harmony,
                stat.by = c(
                  "CNV_score"
                ),
                legend.position = "right", legend.direction = "vertical",
                group.by = "CellType",
                palcolor = zzm60colors[1:9], ylab = "", xlab = "",
                comparisons = list(c("0", "lymphocytes"), c("1", "lymphocytes"),
                                   c("2", "lymphocytes"), c("3", "lymphocytes"),
                                   c("4", "lymphocytes"), c("5", "lymphocytes"),
                                   c("6", "lymphocytes"), c("7", "lymphocytes"))
)
dev.off()

load("./data/sccancer_stemness.Rdata")
sce_harmony$Stemness_score = sccancer_score

pdf("step2_9.pdf")
FeatureStatPlot(sce_harmony,
                stat.by = c(
                  "Stemness_score"
                ),
                legend.position = "right", legend.direction = "vertical",
                group.by = "CellType",
                palcolor = zzm60colors[1:9], ylab = "", xlab = "",
                comparisons = list(c("0", "lymphocytes"), c("1", "lymphocytes"),
                                   c("2", "lymphocytes"), c("3", "lymphocytes"),
                                   c("4", "lymphocytes"), c("5", "lymphocytes"),
                                   c("6", "lymphocytes"), c("7", "lymphocytes"))
)
dev.off()

#################################################################

### Nucleocytoplasmic homologous genes，NCG
library(Seurat)
epi = readRDS("step2data.rds")

epi$Cancer_cells = ifelse(epi$seurat_clusters %in% c(0,3,4,5,7), "cancer", "nc")

epi_sn <- subset(epi, Technology == "snRNA-seq")
epi_sc <- subset(epi, Technology == "scRNA-seq")

diff_DEG_sc <- FindMarkers(epi_sc, min.pct = 0.25, slot = "data", logfc.threshold = 0.25,
                        group.by = "Cancer_cells", ident.1 = "cancer", ident.2 = "nc")

diff_DEG_sn <- FindMarkers(epi_sn, min.pct = 0.25, slot = "data", logfc.threshold = 0.25,
                           group.by = "Cancer_cells", ident.1 = "cancer", ident.2 = "nc")

save(diff_DEG_sc, diff_DEG_sn, file = "sn_sc_degs.Rdata")

load("data/sn_sc_degs.Rdata")
scf1 = diff_DEG_sc %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25) %>% 
  filter(abs(avg_log2FC) > 0.25)

sc_sn_genes = diff_DEG_sn %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25) %>% 
  filter(abs(avg_log2FC) > 0.25) %>% 
  filter(symbol %in% scf1$rowname) %>% 
  dplyr::select(symbol)

a = diff_DEG_sc %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25)

b = diff_DEG_sn %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25)

list_all = list(`scRNA-seq DEGs` = a$symbol, `snRNA-seq DEGs` = b$symbol)
myCol <- c("#0000CC", "#CC0000")

pdf(file = "./Figure/step7_9.pdf", width = 8)
list_all %>% 
  ggvenn(show_elements = F,label_sep = ",",
         digits = 1,
         set_name_color = myCol,
         fill_color = myCol,
         stroke_color = "white",
         set_name_size = 8,
         text_size = 8,
         show_percentage = T
  )
dev.off()

a <- 315
b <- 1094
inter <- 46
phyper(45, b, 50000-b, a, lower.tail = F)
