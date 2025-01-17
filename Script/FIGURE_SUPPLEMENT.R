### Figure supplement

library(ggsci)


zzm9colors <- c('#d4de9c','#94c58f','#86c7b4','#9cd2ed','#a992c0',
                '#ea9994','#f2c396','#bb82b1','#2d3462')

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




col <- pal_aaas()(9)

pdf(file="./Figure/step1_3.pdf", width=15)
CellDimPlot(                                                                           
 srt = sce_harmony, group.by = c("Histopathology", "Technology", "Group", "CellType"),
 reduction = "tsne", palcolor = col                    
)
dev.off()

library(viridis)
library(ggpointdensity)
library(ggplot2)
library(Nebulosa)

sce_N = subset(sce_harmony, Group == "NC")
data = cbind(Embeddings(object=sce_N[['tsne']]),FetchData(sce_N,'Group'))

pdf("step1_4.pdf")
tiff(file = "step1_4.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
ggplot(data = data, 
       mapping = aes(x = tSNE_1, y = tSNE_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()

sce_SCC = subset(sce_harmony, Group == "SCC")
data = cbind(Embeddings(object=sce_SCC[['tsne']]),FetchData(sce_SCC,'Group'))

pdf("step1_5.pdf")
tiff(file = "step1_5.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
ggplot(data = data, 
       mapping = aes(x = tSNE_1, y = tSNE_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()

sce_SCC_1 = subset(sce_harmony, Group == "CCI")
data = cbind(Embeddings(object=sce_SCC_1[['tsne']]),FetchData(sce_SCC_1,'Group'))

pdf("step1_6.pdf")
tiff(file = "step1_6.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
ggplot(data = data, 
       mapping = aes(x = tSNE_1, y = tSNE_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()

sce_SCC_2 = subset(sce_harmony, Group == "CCII")
data = cbind(Embeddings(object=sce_SCC_2[['tsne']]),FetchData(sce_SCC_2,'Group'))

pdf("step1_7.pdf")
tiff(file = "step1_7.tiff", 
     width = 10 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
ggplot(data = data, 
       mapping = aes(x = tSNE_1, y = tSNE_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()

sce_ADC = subset(sce_harmony, Group == "ADC")
data = cbind(Embeddings(object=sce_ADC[['tsne']]),FetchData(sce_ADC,'Group'))

pdf("step1_8.pdf")
tiff(file = "step1_8.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
ggplot(data = data, 
       mapping = aes(x = tSNE_1, y = tSNE_2))+
  NoLegend()+
  geom_pointdensity()+
  scale_color_viridis()+ theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()