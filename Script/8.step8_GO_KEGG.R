rm(list = ls())

##NCG GO/KEGG
library(tidyverse)
library("org.Hs.eg.db")
library(xlsx)
library(RColorBrewer)
library("ggsci")
library("clusterProfiler")
library(ggnewscale)
library(ggfun)
library(grid)
library(ggh4x)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, 
                                     keys=sc_sn_genes$symbol,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="all",
                  readable =T)

go <- as.data.frame(go_kk)

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

kegg = kegg_kk@result
kegg_sub = kegg[1:6,]

go <- go[1:6,] %>% 
  as_tibble() %>% 
  mutate(Count = -1*Count)

go <- go[-1]

data <- rbind(go, kegg_sub)

color = rep("#fec79e", nrow(data))
color[which(data$Count < 0)] <- "#8ec4cb"
data$color <- color

#pdf(file = "./Figure/step2_11.pdf", onefile = F, width = 12)
data %>% 
  ggplot() + 
  geom_col(aes(reorder(Description, Count), Count), fill=color) +
  theme_classic() +
  ylim(-30,30) +
  coord_flip() +
  geom_segment(aes(y=0,yend=0,x=0,xend=13.8)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,size = 25,face = "bold"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=18),
    axis.title.x = element_text(size = 20)
  ) +
  geom_text(data = data[which(data$Count > 0),], aes(x=Description, y=0, label=Description),
            hjust = 1.01, size = 4) +
  geom_text(data = data[which(data$Count < 0),], aes(x=Description, y=0, label=Description),
            hjust = -0.003, size = 4) +
  ggtitle("GO/KEGG Analysis \n FDR < 0.05") +
  scale_x_discrete(expand = (add=c(0,1.5))) +
  geom_segment(aes(y=-1,yend=-20,x=13.5,xend=13.5),
               arrow = arrow(length = unit(0.5,"cm"),type="closed"),
               size=1) +
  geom_segment(aes(y= 1,yend= 20,x=13.5,xend=13.5),
               arrow = arrow(length = unit(0.5,"cm"),type="closed"),
               size=1) +
  annotate("text", x=13.5, y=-22, label="GO",size = 8) +
  annotate("text", x=13.5, y=23, label="KEGG", size = 8)
#dev.off()

### H_vs_L deg GO KEGG
load("./data/HvsL_DEG.Rdata")
DEG = DEG %>% 
  filter(group %in% c("UP", "DOWN"))

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, 
                                     keys=DEG$symbol,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="all",
                  readable =T)

go <- as.data.frame(go_kk)

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

kegg = kegg_kk@result
kegg_sub = kegg[1:6,]

go <- go[1:6,] %>% 
  as_tibble() %>% 
  mutate(Count = -1*Count)

go <- go[-1]

data <- rbind(go, kegg_sub)

color = rep("#fec79e", nrow(data))
color[which(data$Count < 0)] <- "#8ec4cb"
data$color <- color

pdf(file = "./Figure/step3_3.pdf", onefile = F, width = 12)
data %>% 
  ggplot() + 
  geom_col(aes(reorder(Description, Count), Count), fill=color) +
  theme_classic() +
  ylim(-30,30) +
  coord_flip() +
  geom_segment(aes(y=0,yend=0,x=0,xend=13.8)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,size = 25,face = "bold"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=18),
    axis.title.x = element_text(size = 20)
  ) +
  geom_text(data = data[which(data$Count > 0),], aes(x=Description, y=0, label=Description),
            hjust = 1.01, size = 4) +
  geom_text(data = data[which(data$Count < 0),], aes(x=Description, y=0, label=Description),
            hjust = -0.003, size = 4) +
  ggtitle("GO/KEGG Analysis \n FDR < 0.05") +
  scale_x_discrete(expand = (add=c(0,1.5))) +
  geom_segment(aes(y=-1,yend=-20,x=13.5,xend=13.5),
               arrow = arrow(length = unit(0.5,"cm"),type="closed"),
               size=1) +
  geom_segment(aes(y= 1,yend= 20,x=13.5,xend=13.5),
               arrow = arrow(length = unit(0.5,"cm"),type="closed"),
               size=1) +
  annotate("text", x=13.5, y=-22, label="GO",size = 8) +
  annotate("text", x=13.5, y=23, label="KEGG", size = 8)
dev.off()

### new type

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, 
                                     keys=sc_sn_genes$symbol,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

#top 8
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:8)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:8) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

pdf(file = "./Figure/step2_11.pdf", onefile = F, width = 12, height = 10)
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG p.adjust", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP p.adjust", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:16,
                                   labels = plot_df$Description)) +
  ggtitle(label = "GO and KEGG annotation") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(3, "native"),
                    xmax = unit(15, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(42.25, "native"))
dev.off()

###
rm(list = ls())

load("./data/HPV_HvsL_DEG.Rdata")
load("./data/HvsL_DEG.Rdata")

DEG = DEG %>% 
  filter(pvalue < 0.05) %>% 
  filter(log2FoldChange < 0)
DEG = DEG %>% 
  filter(pvalue < 0.05) %>% 
  filter(log2FoldChange > 0)
table(DEG$group)
DEG = DEG %>% 
  filter(group == "DOWN")

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, 
                                     keys=DEG$symbol, #DEG$symbol
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )


go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.9, 
                  qvalueCutoff = 0.9,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

#挑选top 8
GO_top8 <- deg_go_2 %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice(1:8) %>%
  dplyr::ungroup()

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice(1:8) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

pdf(file = "./Figure/step3_6.pdf", onefile = F, width = 12, height = 10)
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#B2DFDB", high = "#4DB6AC", name = "KEGG p.adjust", guide = guide_colorbar(order = 4)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFAB91", high = "#FF7043", name = "MF p.adjust", guide = guide_colorbar(order = 3)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#D1C4E9", high = "#9575CD", name = "CC p.adjust", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFF59D", high = "#FDD835", name = "BP p.adjust", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:32,
                                   labels = plot_df$Description)) +
  ggtitle(label = "GO and KEGG annotation") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 5)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#4DB6AC", "#FF7043", "#9575CD", "#FDD835")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#4DB6AC", "#FF7043", "#9575CD", "#FDD835")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    axis.text.y = element_text(color = rep(c("#4DB6AC", "#FF7043", "#9575CD", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(3, "native"),
                    xmax = unit(15, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(42.25, "native"))
dev.off()

###
geneSet_onco <- read.gmt("./data/c5.all.v2024.1.Hs.symbols.gmt")
geneSet_onco <- read.gmt("./data/c2.all.v2024.1.Hs.symbols.gmt")
geneList <- DEG$log2FoldChange                
names(geneList) <- DEG$symbol     
geneList <- sort(geneList, decreasing = T)  

GSEA_enrichment <- GSEA(geneList,                
                        TERM2GENE = geneSet_onco, 
                        pvalueCutoff = 0.9,      
                        minGSSize = 20,           
                        maxGSSize = 1000,        
                        eps = 0,                 
                        pAdjustMethod = "BH")     

result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)

library(enrichplot)

pdf(file = "./Figure/step4_4_up.pdf", onefile = F, width = 10)
gseaplot2(GSEA_enrichment, "GOBP_HUMORAL_IMMUNE_RESPONSE", color = "red3", pvalue_table = T)
dev.off()

pdf(file = "./Figure/step4_5_down.pdf", onefile = F, width = 20, height = 15)
gseaplot2(GSEA_enrichment, c("GOBP_HUMORAL_IMMUNE_RESPONSE", 
                             "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
                             "GOBP_IMMUNE_EFFECTOR_PROCESS",
                             "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS"), 
          base_size = 20,
          subplots = 1:3,
          color = "red3", pvalue_table = T)
dev.off()

library(GseaVis)

terms = c("GOBP_HUMORAL_IMMUNE_RESPONSE",
          "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
          "GOBP_IMMUNE_EFFECTOR_PROCESS",
          "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
          "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
          "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS")

lapply(terms, function(x){
  gseaNb(object = GSEA_enrichment,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.65,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> gseaList1

#combine
pdf(file = "./Figure/step4_8.pdf", onefile = F, width = 10)
cowplot::plot_grid(plotlist = gseaList1,ncol = 2,align = 'hv')
dev.off()

