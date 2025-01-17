## cell proportion

rm(list = ls())

library(Seurat)
library(tidyverse)
library(speckle)
library(CellBench)
library(limma)
library(ggrepel)
library(ggsci)

memory.limit(size = size_in_megabytes * 1024^2)

sce_harmony = readRDS("./data/step1data.rds")

sce_N_SCC = subset(sce_harmony, Group == "NC"| Group == "SCC")

props <- getTransformedProps(sce_N_SCC$CellType, sce_N_SCC$Histopathology, transform="logit")
group <- c(rep("SCC", 3), rep("NC", 3))
design <- model.matrix(~ 0 + group)
design
mycontr <- makeContrasts(groupSCC-groupNC, levels=design)
#propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, sort=TRUE)
props = as.data.frame(props[3]$Proportions) %>% 
  pivot_wider(names_from="sample", values_from = "Freq") %>%
  as.data.frame()
rownames(props) <- props$clusters
props <- props[-1]
df.fit <- lmFit(props, design)
fit <- contrasts.fit(df.fit, mycontr)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
saveRDS(tempOutput, file="./data/SCC_vs_NC.rds")

#############################################################
tempOutput <- readRDS("./Rdata/SCC_vs_NC.rds")

col <- pal_aaas()(6)
tempOutput <- tempOutput %>% 
  rownames_to_column()
tempOutput$rowname <- factor(tempOutput$rowname, levels = c("lymphocytes",
                                                            "Epithelial cells",
                                                            "Myeloid cells",
                                                            "Fibroblasts",
                                                            "Endothelial cells",
                                                            "Smooth muscle cells"))
rownames(tempOutput) <- tempOutput$rowname
tempOutput <- tempOutput[c("lymphocytes",
                           "Epithelial cells",
                           "Myeloid cells",
                           "Fibroblasts",
                           "Endothelial cells",
                           "Smooth muscle cells"),]

tempOutput_1 = tempOutput %>% 
  dplyr::filter(logFC < 0)

rownames(tempOutput_1) <- tempOutput_1$rowname

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


pdf(file = "./Figure/step1_9.pdf", width = 10, height = 8)
tempOutput %>% 
    ggplot(aes(logFC, -log2(P.Value))) +
      geom_point(aes(size = -log2(P.Value), color = rowname)) +
      theme_classic() +
      theme(legend.position = 'none') +
      theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) +
      geom_vline(xintercept=0, linetype=8, size = 2, color = "grey") +
      geom_hline(yintercept=0, linetype=8, size = 2, color = "grey") +
      scale_x_continuous(name = "log2 fold change")+
      scale_y_continuous(name = "Cluster abundence\n-log P value") + 
      plot.format +
      scale_color_manual(values = col) +
  annotate(geom = "text", x = 0.05, y = 0.2, label = "SCC vs NC Samples", size = 6) +
  geom_text_repel(data = tempOutput_1,
                  aes(x = logFC, y = -log2(P.Value), label = rowname),
                  box.padding = 0.3,
                  nudge_x = 0.02,
                  nudge_y = 0.05,
                  segment.curvature = -0.05,
                  segment.ncp = 3,
                  #segment.angle = 10,
                  direction = "y", 
                  hjust = 0, 
                  color = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
                  size = 5
  ) + 
  geom_text_repel(data =  tempOutput %>% 
                    dplyr::filter(logFC > 0),
                  aes(x = logFC, y = -log2(P.Value), label = rowname),
                  box.padding = 0.3,
                  nudge_x = 0.02,
                  nudge_y = 0.1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  #segment.angle = 20,
                  direction = "y", 
                  hjust = 1,
                  size = 5, 
                  color = c("#631879FF", "#008280FF", "#BB0021FF")
  )

dev.off()
################################################################
sce_CCI_CCII = subset(sce_harmony, Group == "CCI"| Group == "CCII")

props <- getTransformedProps(sce_CCI_CCII$CellType, sce_CCI_CCII$Histopathology, transform="logit")
group <- c(rep("CCII", 3), rep("CCI", 3))
design <- model.matrix(~ 0 + group)
design
mycontr <- makeContrasts(groupCCII-groupCCI, levels=design)
#propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, sort=TRUE)
props = as.data.frame(props[3]$Proportions) %>% 
  pivot_wider(names_from="sample", values_from = "Freq") %>%
  as.data.frame()
rownames(props) <- props$clusters
props <- props[-1]
df.fit <- lmFit(props, design)
fit <- contrasts.fit(df.fit, mycontr)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
saveRDS(tempOutput, file="./data/CCII_vs_CCI.rds")

---------------------------------------------------------------------------------------------
rm(list = ls())
tempOutput <- readRDS("./data/CCII_vs_CCI.rds")

col <- pal_aaas()(6)
tempOutput <- tempOutput %>% 
  rownames_to_column()
tempOutput$rowname <- factor(tempOutput$rowname, levels = c("lymphocytes",
                                                            "Epithelial cells",
                                                            "Myeloid cells",
                                                            "Fibroblasts",
                                                            "Endothelial cells",
                                                            "Smooth muscle cells"))
rownames(tempOutput) <- tempOutput$rowname
tempOutput <- tempOutput[c("lymphocytes",
                           "Epithelial cells",
                           "Myeloid cells",
                           "Fibroblasts",
                           "Endothelial cells",
                           "Smooth muscle cells"),]

pdf(file = "./Figure/step1_10.pdf", width = 10, height = 8)
tempOutput %>% 
  ggplot(aes(logFC, -log2(P.Value))) +
  geom_point(aes(size = -log2(P.Value), color = rowname)) +
  theme_classic() +
  theme(legend.position = 'none') +
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) +
  geom_vline(xintercept=0, linetype=8, size = 2, color = "grey") +
  geom_hline(yintercept=0, linetype=8, size = 2, color = "grey") +
  scale_x_continuous(name = "log2 fold change")+
  scale_y_continuous(name = "Cluster abundence\n-log P value") + 
  plot.format +
  scale_color_manual(values = col) +
  annotate(geom = "text", x = 0.05, y = 0.2, label = "CCII vs CCI Samples", size = 6) +
  geom_text_repel(data = tempOutput %>% 
                    filter(logFC < 0) %>% 
                    filter(rowname == "Fibroblasts"),
                  aes(x = logFC, y = -log2(P.Value), label = rowname),
                  box.padding = 0.5,
                  nudge_x = 0.01,
                  nudge_y = -0.1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  #segment.angle = 10,
                  direction = "y", 
                  hjust = 0, 
                  color = c("#631879FF"),
                  size = 5
  ) + geom_text_repel(data = tempOutput %>% 
                      filter(logFC < 0) %>% 
                      filter(rowname %in% c("lymphocytes", "Myeloid cells")),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = -0.07,
                      nudge_y = 0.1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      #segment.angle = 10,
                      direction = "y", 
                      hjust = 0, 
                      color = c("#3B4992FF", "#008B45FF"),
                      size = 5
  ) + geom_text_repel(data =  tempOutput %>% 
                      dplyr::filter(logFC < 0) %>% 
                      filter(rowname %in% c("Endothelial cells", "Smooth muscle cells")),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = 0.06,
                      nudge_y = 0.1,
                      segment.curvature = 0.1,
                      segment.ncp = 3,
                      #segment.angle = 20,
                      direction = "y", 
                      hjust = 1,
                      size = 5, 
                      color = c("#008280FF", "#BB0021FF")
  ) + geom_text_repel(data =  tempOutput %>% 
                        dplyr::filter(logFC > 0) %>% 
                        filter(rowname == "Epithelial cells"),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = -0.02,
                      nudge_y = 0.1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      #segment.angle = 20,
                      direction = "y", 
                      hjust = 0,
                      size = 5, 
                      color = c("#EE0000FF")
  ) 

dev.off()

###############################################################################
sce_ADC_NC = subset(sce_harmony, Group == "ADC"| Group == "NC")

props <- getTransformedProps(sce_ADC_NC$CellType, sce_ADC_NC$Histopathology, transform="logit")
group <- c(rep("ADC", 3), rep("NC", 3))
design <- model.matrix(~ 0 + group)
design
mycontr <- makeContrasts(groupADC-groupNC, levels=design)
#propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, sort=TRUE)
props = as.data.frame(props[3]$Proportions) %>% 
  pivot_wider(names_from="sample", values_from = "Freq") %>%
  as.data.frame()
rownames(props) <- props$clusters
props <- props[-1]
df.fit <- lmFit(props, design)
fit <- contrasts.fit(df.fit, mycontr)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
saveRDS(tempOutput, file="./data/ADC_vs_NC.rds")
---------------------------------------------------------------------------

rm(list = ls())
tempOutput <- readRDS("./data/ADC_vs_NC.rds")

col <- pal_aaas()(6)
tempOutput <- tempOutput %>% 
  rownames_to_column()
tempOutput$rowname <- factor(tempOutput$rowname, levels = c("lymphocytes",
                                                            "Epithelial cells",
                                                            "Myeloid cells",
                                                            "Fibroblasts",
                                                            "Endothelial cells",
                                                            "Smooth muscle cells"))
rownames(tempOutput) <- tempOutput$rowname
tempOutput <- tempOutput[c("lymphocytes",
                           "Epithelial cells",
                           "Myeloid cells",
                           "Fibroblasts",
                           "Endothelial cells",
                           "Smooth muscle cells"),]

pdf(file = "./Figure/step1_11.pdf", width = 10, height = 8)
tempOutput %>% 
  ggplot(aes(logFC, -log2(P.Value))) +
  geom_point(aes(size = -log2(P.Value), color = rowname)) +
  theme_classic() +
  theme(legend.position = 'none') +
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) +
  geom_vline(xintercept=0, linetype=8, size = 2, color = "grey") +
  geom_hline(yintercept=0, linetype=8, size = 2, color = "grey") +
  scale_x_continuous(name = "log2 fold change")+
  scale_y_continuous(name = "Cluster abundence\n-log P value") + 
  plot.format +
  scale_color_manual(values = col) +
  annotate(geom = "text", x = 0.15, y = 0.5, label = "ADC vs NC Samples", size = 6) +
  geom_text_repel(data = tempOutput %>% 
                    filter(logFC < 0) %>% 
                    filter(rowname == "lymphocytes"),
                  aes(x = logFC, y = -log2(P.Value), label = rowname),
                  box.padding = 0.5,
                  nudge_x = 0.01,
                  nudge_y = -0.1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  #segment.angle = 10,
                  direction = "y", 
                  hjust = 0, 
                  color = c("#3B4992FF"),
                  size = 5
  ) + geom_text_repel(data = tempOutput %>% 
                        filter(logFC < 0) %>% 
                        filter(rowname %in% c("Fibroblasts", "Endothelial cells", "Smooth muscle cells")),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = -0.07,
                      nudge_y = 0.1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      #segment.angle = 10,
                      direction = "y", 
                      hjust = 0, 
                      color = c("#631879FF", "#008280FF", "#BB0021FF"),
                      size = 5
  ) + geom_text_repel(data =  tempOutput %>% 
                        dplyr::filter(logFC > 0) %>% 
                        filter(rowname %in% c("Myeloid cells")),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = 0.06,
                      nudge_y = 0.1,
                      segment.curvature = 0.1,
                      segment.ncp = 3,
                      #segment.angle = 20,
                      direction = "y", 
                      hjust = 1,
                      size = 5, 
                      color = c("#008B45FF")
  ) + geom_text_repel(data =  tempOutput %>% 
                        dplyr::filter(logFC > 0) %>% 
                        filter(rowname == "Epithelial cells"),
                      aes(x = logFC, y = -log2(P.Value), label = rowname),
                      box.padding = 0.5,
                      nudge_x = -0.02,
                      nudge_y = -0.3,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      #segment.angle = 20,
                      direction = "y", 
                      hjust = 0,
                      size = 5, 
                      color = c("#EE0000FF")
  ) 

dev.off()


#################################################################
Cellratio = prop.table(table(sce_harmony$CellType, sce_harmony$Histopathology), margin = 2)
Cellratio <- as.data.frame(Cellratio)
saveRDS(Cellratio, file="./data/allsample_cellpro.rds")

Cellratio = readRDS("./Rdata/allsample_cellpro.rds")
Cellratio = as.data.frame(Cellratio)
colnames(Cellratio) = c("CellType", "Sample", "Freq")

Cellratio$Sample <- factor(Cellratio$Sample, levels = c('ADC1', 'ADC2', 'ADC3',
                                                        'NC1', 'NC2', 'NC3',
                                                        'SCC1', 'SCC2', 'SCC3',
                                                        'CCI1', 'CCI2', 'CCI3',
                                                        'CCII1', 'CCII2', 'CCII3'))
Cellratio$CellType = factor(Cellratio$CellType, levels = c("lymphocytes",
                                                           "Epithelial cells",
                                                           "Myeloid cells",
                                                           "Fibroblasts",
                                                           "Endothelial cells",
                                                           "Smooth muscle cells") )


pdf(file = "./Figure/step1_12.pdf", width = 10, height = 8)
ggplot(Cellratio,aes(Sample,Freq,fill = CellType)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  #coord_flip() +
  theme(legend.position = 'none') +
  plot.format + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = col)
dev.off()

###########################################################################
## cor

rm(list = ls())

library(Seurat)
library(ggplot2)
library(tidyverse)
library(BiocParallel)
library(AnnoProbe)
library(SCP)

register(MulticoreParam(workers = 8, progressbar = TRUE))
memory.limit(size = size_in_megabytes * 1024^2)

sce_harmony = readRDS("step1data.rds")

pdf("./Figure/step1_13.pdf", width=12, height=10)
ht <- CellCorHeatmap(
  srt_query = sce_harmony,
  query_group = "CellType",
  #query_assay = "data",
  #query_reduction = "harmony",
  query_dims = 12,
  #exp_legend_title = "",
  distance_metric = "pearson",
  nlabel = 3, 
  label_by = "row",
  row_title = "",
  column_title = "",
  heatmap_palcolor =  c("#FFF0F5","#FFC0CB","#DC143C"),
  query_group_palcolor = c("#3B4992FF","#EE0000FF","#008B45FF","#631879FF","#008280FF","#BB0021FF","#5F559BFF", "#A20056FF","#808180FF"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  label_size = 15,
  units="mm"
)
ht$plot
dev.off()

