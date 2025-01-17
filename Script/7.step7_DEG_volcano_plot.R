### TCGA DEGS
rm(list = ls())
library(ggrepel)
library(patchwork)
library(tidyverse)

load("data/step5_all_survial_data.Rdata")

cg = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1")
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
group.gsva <- ifelse(dat.gsva < 0, "NCG_low", "NCG_high")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)

library(DESeq2)

raw = all_data[-(1:4)]
raw = as.data.frame(t(raw))

condition <- factor(life_data$group, levels = c("NCG_high","NCG_low"))
colData <- data.frame(row.names=colnames(raw), condition)

dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, 
               contrast=c("condition","NCG_high","NCG_low"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG = na.omit(DEG)
DEG = DEG %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  mutate(group = case_when(padj < 0.05 & log2FoldChange > 1 ~ "UP",
                           padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
                           TRUE ~ 'no significance'))
save(DEG, file="./data/HvsL_DEG.Rdata")
save(DEG, file="./data/HPV_HvsL_DEG.Rdata")

#####

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)
library(ggrastr)

load("./data/HvsL_DEG.Rdata")
a_DEG = DEG
load("./data/HPV_HvsL_DEG.Rdata")
hpv_DEG = DEG
rm(DEG)

a_DEG$log_pval <- -log10(a_DEG$padj)
hpv_DEG$log_pval <- log10(hpv_DEG$padj)

a = a_DEG %>% 
  filter(group != "no significance") %>% 
  filter(symbol %in% markgene2) %>% 
  filter(symbol %in% c('MUC5AC', 'BPIFB1', 'MUC13', 'CXCL10', 'AGR2', 'EGF',
                       'IDO2',  'TFF3', 'PIGR', 'CCL19', 'MMP7', 'HHLA2', 'FOSB',
                       'MAGEA4', 'KRT1', 'CFTR', 'MUC5B',  'MS4A8', 'AGR3',
                       'TSPAN8', 'PROM1', 'NTS', 'IRX4'))

b = hpv_DEG %>% 
  filter(group != "no significance") %>% 
  filter(symbol %in% markgene2) %>% 
  filter(symbol %in% c("LAG3", "CXCL10", "LTB", "EGF", "CXCL11", "NKG7",
                       "GZMH", "CCL19", "ISG15", "CCL5", "PDCD1", "GZMK",
                       "GNLY", "HLA-DQB2", "CD200R1",  "HLA-DQA2", "MAGEA4",
                       "GZMB", "CXCL9", "IFNG", "HHLA2", "TSPAN8", "MS4A8",
                       "CRNN", "TFF3", "SBSN", "BPIFB1"))

df <- rbind(a_DEG,hpv_DEG)

rect1 <- rasterGrob(matrix(rev(alpha(c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1','#DD3497', '#7A0177'),0.5)),nrow = 1),                    width=unit(1,"npc"), height = unit(1,"npc"), interpolate = TRUE)
rect2<- rasterGrob(matrix(alpha(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d'),0.5),nrow = 1),                    width=unit(1,"npc"), height = unit(1,"npc"), interpolate = TRUE)
rect3<- rasterGrob(matrix(rev(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4")),nrow = 1), width=unit(1,"npc"),
                   height = unit(1,"npc"), interpolate = TRUE) 
rect4<- rasterGrob(matrix(rev(c("#DC6F58","#F7B698","#FAE7DC","#E1EDF3")),nrow = 1),width=unit(1,"npc"),
                   height = unit(1,"npc"), interpolate = TRUE)

p = ggplot(df, aes(log2FoldChange,log_pval))+
  annotation_custom(rect1,xmin=-Inf,xmax=-1,ymin=1.3,ymax=Inf)+
  annotation_custom(rect2, xmin=1,xmax=Inf,ymin=1.3,ymax=Inf)+
  annotation_custom(rect3,xmin=-Inf,xmax=-1,ymin=-1.3,ymax=-Inf)+
  annotation_custom(rect4, xmin=1,xmax=Inf,ymin=-1.3,ymax=-Inf)+
  geom_hline(aes(yintercept=1.3), color = "#999999", linetype="dashed", size=1) +
  geom_hline(aes(yintercept=-1.3), color = "#999999", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=-1), color = "#999999", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=1), color = "#999999", linetype="dashed", size=1) + 
  geom_point_rast(stroke = 0.5, size=2, shape=21, fill="grey") +
  geom_point_rast(data = df[abs(df$log_pval) >= 1.3  &  abs(df$log2FoldChange) >=1,],
                  aes(fill=log2FoldChange),stroke = 0.5, size=3, shape=21)+ 
  scale_fill_gradient2(limits=c(-max(df$log2FoldChange), max(df$log2FoldChange)),
                       low="blue", mid="whitesmoke", high = "red", na.value = 'blue')+ 
  geom_text_repel(data = a[a$log2FoldChange >0,], aes(label=symbol), 
                  nudge_y = 25, 
                  nudge_x = 3, 
                  color = 'black',
                  size = 4,fontface = 'italic',
                  min.segment.length = 0,
                  max.overlaps = Inf, 
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm')
                  ) +
  geom_text_repel(data = a[a$log2FoldChange <0,], aes(label=symbol), 
                  nudge_y = 25, 
                  nudge_x = -3, 
                  color = 'black',
                  size = 4,fontface = 'italic',
                  min.segment.length = 0,
                  max.overlaps = Inf, 
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm')
                  ) +
  geom_text_repel(data = b[b$log2FoldChange >0,], aes(label=symbol), 
                  nudge_y = -25, 
                  nudge_x = 3, 
                  color = 'black',
                  size = 4,fontface = 'italic',
                  min.segment.length = 0,
                  max.overlaps = Inf, 
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm')
  ) +
  geom_text_repel(data = b[b$log2FoldChange <0,], aes(label=symbol), 
                  nudge_y = -25, 
                  nudge_x = -3, 
                  color = 'black',
                  size = 4,fontface = 'italic',
                  min.segment.length = 0,
                  max.overlaps = Inf, 
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm')
  )+
  theme_classic()+ 
  theme(aspect.ratio = 1,
        axis.text.x=element_text(colour="black",size =12),
        axis.text.y=element_text(colour="black",size =12),
        axis.ticks=element_line(colour="black"),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,1,0),"line"), 
        legend.direction = 'horizontal', 
        legend.position = 'top',
        legend.justification=c(1,2),
        legend.key.width=unit(0.5,"cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.title.position = 'top')+
  geom_rect(aes(xmin =-Inf, xmax = Inf, ymin = 1.3, ymax = 40),
            fill = "transparent", color = "red", size = 1.5,linetype="dashed")+
  geom_rect(aes(xmin =-Inf, xmax = Inf, ymin = -1.3, ymax = -40),
            fill = "transparent", color = "blue", size = 1.5,linetype="dashed")+
  scale_y_continuous(limits = c(-40, 40), expand = c(0,0),breaks = c(-40,0,40), labels = c(40,0,40))+
  scale_x_continuous(limits = c(-10, 10))


pdf(file = "./Figure/step6_3.pdf", onefile = F, width = 8, height = 7.7)
ggdraw(xlim = c(0, 1), ylim = c(0,1.1)) +
  draw_plot(p, x = 0, y =0)+
  draw_line(x = c(0.3,0.8), y = c(0.01,0.01),lineend = "round", 
            size =1, col = "black", 
            arrow=arrow(angle = 15, length = unit(0.1,"inches"),type = "closed"))+
  draw_text(text = "Log2FC", size = 12,
            x = 0.2, y = 0.02,color="black",fontface = "italic")+ 
  draw_line(x = c(0.1,0.1), y = c(0.2,0.85),
            lineend = "round",size =1, col = "black",   
            arrow=arrow(angle = 15, length = unit(0.1,"inches"),type = "closed"))+
  draw_text(text = "-Log(P)", size = 12,angle=90,
            x = 0.1, y = 0.15,color="black",fontface = "italic")+ 
  draw_text(text = "HRG vs LRG", size = 12,angle=90,x = 0.9, y = 0.7,color="black",fontface = "bold")+
  draw_text(text = "HPV+ HRG vs HPV+ LRG", size = 12,angle=90,x = 0.9, y = 0.3,color="black",fontface = "bold")
dev.off()


#volcano
load("./data/HvsL_DEG.Rdata")

DEG$group <- factor(DEG$group, levels = c("UP","DOWN","no significance"))

mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggplot(data = DEG,
       aes(x = log2FoldChange, y = -log10(padj), color = group)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-11, 11), breaks = seq(-12, 12, by = 4))+
  scale_y_continuous(expand = expansion(add = c(1, 0)),
                     limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  mytheme

DEG$curve_y <- case_when(
  DEG$log2FoldChange > 0 ~ 1/(DEG$log2FoldChange-1) + (-log10(0.05)),
  DEG$log2FoldChange <= 0 ~ 1/(-DEG$log2FoldChange-1) + (-log10(0.05))
)

DEG$group2 <- case_when(
  -log10(DEG$padj) > DEG$curve_y & DEG$log2FoldChange >= 1 ~ 'up',
  -log10(DEG$padj) > DEG$curve_y & DEG$log2FoldChange <= -1 ~ 'down',
  TRUE ~ 'no significance'
)

DEG$group2 <- factor(DEG$group2, levels = c("up","down","no significance"))

f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(0.05))
  dff <- rbind(data.frame(x = inputx + 1, y = y),
               data.frame(x = -(inputx + 1), y = y))
  return(dff)
}

dff_curve <- f(5)
head(dff_curve)

ggplot(data = DEG,
       aes(x = log2FoldChange, y = -log10(padj), color = group2)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-11, 11), breaks = seq(-12, 12, by = 4))+
  scale_y_continuous(expand = expansion(add = c(1, 0)),
                     limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y), 
            color = "black",lty = "dashed", size = 0.7) +
  mytheme


top10 <- filter(DEG, group2 != "none") %>%
  distinct(Symbol, .keep_all = T) %>%
  top_n(10, abs(log2FoldChange))

p5 <- p3 +
  geom_text_repel(data = top10,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)

###
DEG = diff_DEG_sc %>% 
  rownames_to_column("symbol") %>% 
  as_tibble()

DEG$curve_y <- case_when(
  DEG$avg_log2FC > 0 ~ 1/(DEG$avg_log2FC-0.25) + (-log10(0.05)),
  DEG$avg_log2FC <= 0 ~ 1/(-DEG$avg_log2FC-0.25) + (-log10(0.05))
)

DEG$group <- case_when(
  -log10(DEG$p_val_adj) > DEG$curve_y & DEG$avg_log2FC>= 0.25 ~ 'up',
  -log10(DEG$p_val_adj) > DEG$curve_y & DEG$avg_log2FC <= -0.25 ~ 'down',
  TRUE ~ 'no significance'
)

DEG$group <- factor(DEG$group, levels = c("up","down","no significance"))

f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(0.05))
  dff <- rbind(data.frame(x = inputx + 0.25, y = y),
               data.frame(x = -(inputx + 0.25), y = y))
  return(dff)
}

dff_curve <- f(5)
head(dff_curve)

cg = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1")

pdf(file = "./Figure/step2_12.pdf", onefile = F, width = 8)
ggplot(data = DEG,
       aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2))+
  scale_y_continuous(expand = expansion(add = c(1, 0)),
                     limits = c(0, 300), breaks = seq(0, 300, by = 100)) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y), 
            color = "black",lty = "dashed", size = 0.7) +
  mytheme +
  geom_text_repel(data = diff_DEG_sc %>% 
                    rownames_to_column("symbol") %>% 
                    as_tibble() %>% 
                    filter(p_val_adj < 0.05) %>% 
                    filter(pct.1 > 0.25) %>% 
                    filter(pct.2 > 0.25) %>% 
                    filter(symbol %in% cg),
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 1, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.1,
                  nudge_x = 0,
                  nudge_y = 1)
dev.off()

###

set.seed(123)  
log_min <- log10(1e-100) 
log_max <- log10(1e-50) 

log_random_values <- runif(2, min = log_min, max = log_max)

log_random_values_transformed <- log_random_values^1 

random_values <- 10^log_random_values_transformed


DEG_1 = diff_DEG_sn %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj == 0) %>% 
  filter(symbol %in% cg) %>% 
  mutate(p_val_adj = p_val_adj + random_values)

DEG_2 = diff_DEG_sn %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj != 0)

DEG = rbind(DEG_1, DEG_2)

DEG$curve_y <- case_when(
  DEG$avg_log2FC > 0 ~ 1/(DEG$avg_log2FC-0.25) + (-log10(0.05)),
  DEG$avg_log2FC <= 0 ~ 1/(-DEG$avg_log2FC-0.25) + (-log10(0.05))
)

DEG$group <- case_when(
  -log10(DEG$p_val_adj) > DEG$curve_y & DEG$avg_log2FC>= 0.25 ~ 'up',
  -log10(DEG$p_val_adj) > DEG$curve_y & DEG$avg_log2FC <= -0.25 ~ 'down',
  TRUE ~ 'no significance'
)

DEG$group <- factor(DEG$group, levels = c("up","down","no significance"))

f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(0.05))
  dff <- rbind(data.frame(x = inputx + 0.25, y = y),
               data.frame(x = -(inputx + 0.25), y = y))
  return(dff)
}

dff_curve <- f(5)
head(dff_curve)

cg = c("ZFP36L1", "P4HB", "TTC3", "TPM3", "TPT1", "TM9SF2", "ACTR3", "MGST1")

pdf(file = "./Figure/step2_13.pdf", onefile = F, width = 8)
ggplot(data = DEG,
       aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1))+
  scale_y_continuous(expand = expansion(add = c(1, 0)),
                     limits = c(0, 300), breaks = seq(0, 350, by = 100)) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y),
            color = "black",lty = "dashed", size = 0.7) +
  mytheme +
  geom_text_repel(data = DEG %>% 
                    filter(p_val_adj < 0.05) %>% 
                    filter(pct.1 > 0.25) %>% 
                    filter(pct.2 > 0.25) %>% 
                    filter(symbol %in% cg),
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 1, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.1,
                  nudge_x = 0,
                  nudge_y = 1)
dev.off()


#####
rm(list = ls())
load("data/step2.Rdata")

Tumor_count_data = as.data.frame(Tumor_count_data)
Normal_count_data = as.data.frame(Normal_count_data)

dat = cbind(Normal_count_data, Tumor_count_data)

nor_data = dat %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble() %>% 
  inner_join(gene_id_pc, by = "gene_id") %>% 
  as.data.frame()

nb = length(dat) + 1

nor_data$mean_expression <- rowMeans(nor_data[2:nb])
nor_data <- nor_data[order(nor_data$gene_name, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$gene_name), ]
rownames(nor_data) = nor_data$gene_name
nor_data = nor_data[2:nb]

all_data = as.data.frame(t(nor_data))

data = all_data[colnames(all_data) %in% sc_sn_genes$symbol]
data = log2(data + 1)
data$group = c(rep("Normal",3),rep("Cancer", 306))

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


tiff(file = "./Figure/step6_1.tiff", 
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
        axis.title.x = element_blank(),
        ) +
  plot.format +
  scale_fill_manual(values = c("Cancer" = "#DDDDDD", "Normal" = "#E6B745")) +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("TCGA-CESC(n=309)")
dev.off()
