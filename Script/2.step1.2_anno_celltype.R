## annotation celltype
rm(list = ls())

library(Seurat)
library(tidyverse)
library(SCP)
library(BiocParallel)
library(SingleR)
library(celldex)

register(MulticoreParam(workers = 8, progressbar = TRUE))
sce_harmony = readRDS("./data/step1data.rds")

aa <- as.character(sce_harmony$RNA_snn_res.0.1)
names(aa) <- rownames(sce_harmony@meta.data)
sce_harmony@active.ident <- factor(aa)

pdf(file="t1.pdf",width=30,height=11)
clusterMarker=c("EPCAM", "KRT7",
                  "EMCN", "PLVAP",
                  "KIT", "CPA3",
                  "TAGLN", "ACTA2", "MYH11",
                  "PDGFRA", "LAMA2",
                  "CD3E", "CD3D", "CD2",
                  "MZB1", "JCHAIN",
                  "CD19", "MS4A1",
                  "CD14", "C1QA")
DotPlot(object = sce_harmony, features = features)
dev.off()

VlnPlot(epi,
        features = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1"),
        split.by = 'seurat_clusters',stack = TRUE, 
        flip = T) + 
  theme(legend.position='none')+labs(x='',y='')



#refdata <- HumanPrimaryCellAtlasData()
#testdata <- GetAssayData(sce_harmony, slot="data")
#clusters <- sce_harmony$RNA_snn_res.0.1
#cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
#                    method = "cluster", clusters = clusters, 
#                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
#celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
#write.csv(celltype,"./data/celltype_singleR.csv",row.names = F)

current.cluster.ids <- c(0:8)
new.cluster.ids <- c("lymphocytes","Epithelial cells","Epithelial cells",
                     "Fibroblasts", "Myeloid cells", "Epithelial cells", 
                     "Smooth muscle cells", "Endothelial cells", "Epithelial cells")

sce_harmony$CellType <- plyr::mapvalues(x = as.integer(as.character(sce_harmony$RNA_snn_res.0.1)), from = current.cluster.ids, to = new.cluster.ids)

sce_harmony$Cluster = sce_harmony$RNA_snn_res.0.1
sce_harmony$Cluster = factor(sce_harmony$Cluster, levels= c("0", "1", "2", "3",
                                                            "4", "5", "6", "7",
                                                            "8")
)

#sce_harmony$Histopathology = factor(sce_harmony$Histopathology, 
#                                    levels =  c('scRNA_N1', 'scRNA_N2', 'scRNA_HPV_N1',
#                                                'scRNA_HPV_N2', 'scRNA_HPV_HSIL1',
#                                                'scRNA_HPV_HSIL2', 'scRNA_HPV_SCC1',
#                                                'scRNA_HPV_SCC2', 'scRNA_HPV_SCC3',
#                                                'scRNA_HPV_SCC4', 'scRNA_HPV_SCC5',
#                                                'snRNA_HPV_HSIL', 'snRNA_HPV_SCC_CCI1',
#                                                'snRNA_HPV_SCC_CCI2', 'snRNA_HPV_SCC_CCII1',
#                                                'snRNA_HPV_SCC_CCII2', 'snRNA_HPV_SCC_CCII3')
#)

sce_harmony$orig.ident[sce_harmony$orig.ident == "N1"] <- "NC1"
sce_harmony$orig.ident[sce_harmony$orig.ident == "N2"] <- "NC2"
sce_harmony$orig.ident[sce_harmony$orig.ident == "N3"] <- "NC3"

sce_harmony$Histopathology = sce_harmony$orig.ident
sce_harmony$Group = case_when(sce_harmony$orig.ident %in% c("ADC1", "ADC2", "ADC3") ~ "ADC", 
                              sce_harmony$orig.ident %in% c("CCI1", "CCI2", "CCI3") ~ "CCI",
                              sce_harmony$orig.ident %in% c("CCII1", "CCII2", "CCII3") ~ "CCII",
                              sce_harmony$orig.ident %in% c("NC1", "NC2", "NC3") ~ "NC",
                              sce_harmony$orig.ident %in% c("SCC1", "SCC2", "SCC3") ~ "SCC")

pdf(file="./Figure/step1_1.pdf",width=15)
ht <- GroupHeatmap(
  srt = sce_harmony,
  features = c(
    "CD3E", "CD3D", "CD2", # lymphocytes
    "KRT5", "KRT8", # Epithelial cells
    "ITGAX", "LYZ", # Myeloid cells
    "DCN", # Fibroblasts
    "VWF", "CDH5", # Endothelial cells
    "ACTG2", "MYH11" # Smooth muscle cells
  ),                        
  group.by = c("CellType", "Cluster"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase"),        
  cell_annotation_palette = c("Dark2"),
  show_row_names = F, row_names_side = "left",
  add_dot = TRUE, add_reticle = F
)
print(ht$plot)
dev.off()

saveRDS(sce_harmony, "step1data.rds")


library(Seurat)
library(patchwork)

# 
genes <- c(
           "CD3E", "CD3D", "CD2", # lymphocytes
           "KRT5", "KRT8", # Epithelial cells
           "ITGAX", "LYZ", # Myeloid cells
           "DCN", # Fibroblasts
           "VWF", "CDH5", # Endothelial cells
           "ACTG2", "MYH11" # Smooth muscle cells
           )  

#pdf("step1_1-1.pdf", width=10, height=10, onefile = FALSE)
FeatureStatPlot(sce_harmony,
                stat.by = genes,
                fill.by = "group",
                plot_type = "violin", 
                group.by = "CellType",
                bg.by = "CellType",
                stack = TRUE,
                flip = FALSE,
                palcolor = col,
                bg_palcolor = col
)
#dev.off()


# t-SNE
plots <- lapply(genes, function(gene) {
  FeaturePlot(sce_harmony, features = gene, reduction = "tsne" ,pt.size = 1) + 
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

combined_plot <- wrap_plots(plots, ncol = 4)

pdf("./Figure/step1_2.pdf", width = 10)
print(combined_plot)
dev.off()

CellType <- c(
              "lymphocytes", "Epithelial cells", "Myeloid cells",
              "Fibroblasts", "Endothelial cells", "Smooth muscle cells"
              )

comparisons <- list(
  c("ADC", "NC"),
  c("SCC", "NC"),
  c("SCC_II", "SCC_I")
)

results <- lapply(CellType, function(cell_type) {
  lapply(comparisons, function(comp) {
    sce <- subset(sce_harmony, Class %in% comp & CellType == cell_type)
    diff_DEG <- FindMarkers(sce, min.pct = 0.25, slot = "data", logfc.threshold = 0.25,
                            group.by = "Class", ident.1 = comp[1], ident.2 = comp[2])
    diff_DEG$CellType <- cell_type
    diff_DEG$Comparison <- paste(comp[1], "vs", comp[2])
    diff_DEG
  })
})

flattened_results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))

head(flattened_results)

load("./Rdata/three_group_celltype_DEGs.Rdata")
flattened_results %>% 
  filter(Comparison == "HSIL vs NC") %>% 
  filter(CellType == "Epithelial cells")

#########################################################################################
#########################################################################################

gene_sets <- list(
  ssec = c("KRT5", "TP63"),
  gec = c("KRT8", "KRT19"),
)

for (pathway_name in names(gene_sets)) {
  scRNA <- AddModuleScore(
    object = scRNA,
    features = list(gene_sets[[pathway_name]]),
    name = pathway_name
  )
}

pdf("t1.pdf")
FeatureStatPlot(scRNA, stat.by = c("KRT5", "KRT8"), group.by = "RNA_snn_res.0.1")
dev.off()

expression_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(expression_matrix) <- paste0("Gene", 1:100)
colnames(expression_matrix) <- paste0("Sample", 1:10)

gene_set <- list(
  GeneSet1 = c("Gene1", "Gene2", "Gene3"),
  GeneSet2 = c("Gene4", "Gene5", "Gene6", "Gene7")
)

ssgsea_scores <- ssgseaParam(expression_matrix, gene_set)
c = gsva(ssgsea_scores)




