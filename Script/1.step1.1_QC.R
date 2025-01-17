rm(list = ls())
#Sys.setenv(R_MAX_NUM_DLLS=999)
memory.limit(size = size_in_megabytes * 1024^2)
options(stringsAsFactors = F)
gc()

library(Seurat)
library(tidyverse)
library(harmony)
library(DoubletFinder)

source("Rfunction_DoubletFinder.R")

folders=list.files('./',pattern='[1-6]$')
names(folders) = c('ADC1', 'ADC2', 'ADC3',
                   'CCI1', 'CCI2', 'CCI3',
                   'CCII1', 'CCII2', 'CCII3',
                   'N1', 'N2', 'N3',
                   'SCC1', 'SCC2', 'SCC3')

scRNAlist <- list()
for(i in 1:length(folders)){
  counts <- Read10X(data.dir = folders[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA)/log10(scRNAlist[[i]]$nCount_RNA)
}

subscRNAlist1 = scRNAlist[c(1:3, 10:15)]
subscRNAlist2 = scRNAlist[4:9]

### QC
sce_sc <- merge(subscRNAlist1[[1]], 
                y = c(subscRNAlist1[[2]],subscRNAlist1[[3]],subscRNAlist1[[4]],
                      subscRNAlist1[[5]],subscRNAlist1[[6]],
                      subscRNAlist1[[7]],subscRNAlist1[[8]],
                      subscRNAlist1[[9]] 
                ))


sce_sn <- merge(subscRNAlist2[[1]], 
                y = c(subscRNAlist2[[2]],subscRNAlist2[[3]],
                      subscRNAlist2[[4]],subscRNAlist2[[5]],
                      subscRNAlist2[[6]] 
                ))

sce_sc$orig.ident = factor(sce_sc$orig.ident, levels=c("N1", "N2", "N3",
                                                       "SCC1", "SCC2", "SCC3",
                                                       "ADC1", "ADC2", "ADC3"))

pdf("QC-1.pdf", width = 15)
VlnPlot(sce_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 4, group.by="orig.ident")
dev.off()

pdf("QC-1.pdf", width = 12)
VlnPlot(sce_sn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 4, group.by="orig.ident")
dev.off()

### del
subscRNAlist1 <- lapply(X = subscRNAlist1, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 &
                log10GenesPerUMI > 0.7 &
                percent.mt < 10 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})

subscRNAlist2 <- lapply(X = subscRNAlist2, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 &
                log10GenesPerUMI > 0.7 &
                percent.mt < 1 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})

scRNAlist = c(subscRNAlist1, subscRNAlist2)

for (i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], npcs = 30, verbose = FALSE)
  scRNAlist[[i]] <- RunUMAP(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindNeighbors(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindClusters(scRNAlist[[i]], resolution = 0.8)
  
  Doubletrate = ncol(scRNAlist[[i]])*8*1e-6
  sweep.res.list <- paramSweep(scRNAlist[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  sweep.stats[order(sweep.stats$BCreal),]
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- modelHomotypic(scRNAlist[[i]]$seurat_clusters)
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  scRNAlist[[i]] <- doubletFinder(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  scRNAlist[[i]]$doubFind_res = scRNAlist[[i]]@meta.data %>% select(contains('DF.classifications'))
  scRNAlist[[i]]$doubFind_score = scRNAlist[[i]]@meta.data %>% select(contains('pANN'))
  scRNAlist[[i]] = subset(scRNAlist[[i]], doubFind_res == "Singlet")
}

sce <- merge(scRNAlist[[1]], 
             y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],
                   scRNAlist[[5]],scRNAlist[[6]],
                   scRNAlist[[7]],scRNAlist[[8]],
                   scRNAlist[[9]],scRNAlist[[10]],
                   scRNAlist[[11]],scRNAlist[[12]],
                   scRNAlist[[13]],scRNAlist[[14]],
                   scRNAlist[[15]] 
             ), 
             #add.cell.ids = folders, 
             project = "scRNA-snRNA")

sce$Technology = ifelse(sce$orig.ident %in% c("CCI1", "CCI2", "CCI3",
                                              "CCII1", "CCII2", "CCII3"),
                        "snRNA-seq",
                        "scRNA-seq"
)

sce <- NormalizeData(sce)

test.sec <- CellCycleScoring(sce, s.features = Seurat::cc.genes.updated.2019$s.genes,
                             g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
sce$Phase = test.sec$Phase
sce$S.Score = test.sec$S.Score
sce$G2M.Score = test.sec$G2M.Score

sce = FindVariableFeatures(sce, selection.method = "vst",nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score"))
sce <- RunPCA(sce, npcs = 50, verbose = FALSE)

pdf("ElbowPlot1.pdf")
ElbowPlot(sce, ndims=50, reduction="pca")
dev.off()

pc.num = 1:12
sce_harmony <- RunHarmony(sce, group.by.vars = c("Technology", "orig.ident"))
sce_harmony <- RunUMAP(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- RunTSNE(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- FindNeighbors(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- FindClusters(sce_harmony, reduction = "harmony", resolution = 0.1)

saveRDS(sce_harmony, "step1data.rds")

