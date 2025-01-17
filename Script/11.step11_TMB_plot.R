### TMB
rm(list = ls())

library(tidyverse)
library(maftools)
library(TCGAbiolinks)
library(GSVA)
library(GSEABase)

query <- GDCquery(
  project = "TCGA-CESC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-CESC_SNP.Rdata")

load("data/TCGA-CESC_SNP.Rdata")
load("./data/step5_all_survial_data.Rdata")

cg = c("TTC3", "MGST1", 'ACTR3')

da = all_data[-(1:4)]
da = log2(da + 1)
da = as.data.frame(t(da))

cg=cg[cg %in% rownames(da)]
cgl = list(cg=cg)

geneSet <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, cgl, names(cgl)))

da = as.matrix(da)
a = gsvaParam(
  da,
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

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)
#life_data$group = factor(life_data$group, levels = c("NCG_low", "NCG_high"))

rownames(life_data) = life_data$ID1
life_data = life_data[-1]
colnames(life_data) = c("fustat", "futime", "Tumor_Sample_Barcode", "Group")

low_life = life_data %>% 
  filter(Group == "LRG")
high_life = life_data %>% 
  filter(Group == "HRG")

CESC_maf = data
CESC_maf$Barcode = substr(CESC_maf$Tumor_Sample_Barcode, 1, 12)
CESC_maf = CESC_maf %>% 
  filter(Barcode %in% life_data$Tumor_Sample_Barcode)
CESC_maf_l = CESC_maf %>% 
  filter(Barcode %in% low_life$Tumor_Sample_Barcode)
CESC_maf_h = CESC_maf %>% 
  filter(Barcode %in% high_life$Tumor_Sample_Barcode)

mafl = read.maf(maf = CESC_maf_l, isTCGA = T)
mafh = read.maf(maf = CESC_maf_h, isTCGA = T)

pt.vs.rt <- mafCompare(m1 = mafh, m2 = mafl, m1Name = 'HRG',
                       m2Name = 'LRG', minMut = 5)

pdf(file = "./Figure/step4_6.pdf", onefile = F, width = 10)
tiff(file = "./Figure/step4_6.tiff", 
     width = 11 * 300,   
     height = 12 * 300,  
     res = 300, 
     compression = "lzw")
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()

cg = c("OBSCN", "ERBB3", "DAAM1", "DHX16", "UPF1", "FAT1", 
             "ARHGAP36", "C6", "CLIP3", "EZR", "KLHL1", "SLC13A1", 
             "DNAH17", "SACS", "SI", "AKT1", "BCAN", "CAMSAP1", 
             "KDM3B", "NAV1", "OTOF", "PCDHA5", "PCDHGA2", "POLRMT", 
             "ROBO2", "TMEM51", "ZNF611", "MXRA5", "C2orf16", "LTBP4")

data = all_data[-c(1:4)][colnames(all_data[-c(1:4)]) %in% cg]
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

pdf(file = "./Figure/step4_7.pdf", onefile = F, width = 8)
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

## live analysis

library("survival")
library("survminer")
library("clusterProfiler")
library(tidyverse)

head(all_data)
hub_gene = c("SI", "EZR", "DAAM1")

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
    pdf(sprintf("live plot TMB %s.pdf", gene), width = 8, onefile = FALSE)
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

###
maf = read.maf(maf = CESC_maf, clinicalData = life_data, isTCGA = T)

group_colors = list(
  Group = c(
    "HRG" = "#DDDDDD",
    "LRG" = "#E6B745"
  )
)

vc_cols = RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'Splice_Site',
  'In_Frame_Ins'
)

ge = getGO("GO:0006959")
humoral_immune_response = c("NOTCH1", "HLA-A", "NOTCH2", "PTPRC", "NOD2", "SVEP1", "CFHR5", "C3", "DMBT1",
                            "KRT1")

pdf(file = "./Figure/step4_9_humoral_immune_response.pdf", onefile = F, width = 10)
oncoplot(
  maf = maf,
  clinicalFeatures = 'Group',
  sortByAnnotation = TRUE,
  genes = humoral_immune_response,
  annotationColor = group_colors,
  colors = vc_cols,
  titleText = "Humoral Immune Response"
)
dev.off()

ge = getGO("GO:0061844")
ahirm = c("NOD2", "SPINK5", "EVPL", "MUC7", "REG1A", "H2BC8", "H2BC7", "H2BC4", "ELANE", "GAPDH")

pdf(file = "./Figure/step4_10_ahirm.pdf", onefile = F, width = 10)
oncoplot(
  maf = maf,
  clinicalFeatures = 'Group',
  sortByAnnotation = TRUE,
  genes = ahirm,
  annotationColor = group_colors,
  colors = vc_cols,
  titleText = "Antimicrobial Humoral Immune Response Mediated By Antimicrobial Peptide"
)
dev.off()

vc_cols = RColorBrewer::brewer.pal(n = 6, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'Splice_Site'
)

pdf(file = "./Figure/step3_16_low.pdf", onefile = F, width = 10)
oncoplot(
  maf = mafl,
  top = 20,
  #draw_titv = TRUE,
  titleText = "Top 20 Mutated Genes of LRG(n=114)",
  colors = vc_cols
)
dev.off()

pdf(file = "./Figure/step3_17_high.pdf", onefile = F, width = 10)
oncoplot(
  maf = mafh,
  top = 20,
  #draw_titv = TRUE,
  titleText = "Top 20 Mutated Genes of HRG(n=100)",
  colors = vc_cols
)
dev.off()

### TMB
pdf(file = "./Figure/step3_18_low.pdf", onefile = F)
maf.TMB <- tcgaCompare(
  maf = maf, cohortName = 'NCG_low', 
  logscale = TRUE
)
dev.off()

pdf(file = "./Figure/step3_19_high.pdf", onefile = F)
maf.TMB <- tcgaCompare(
  maf = maf, cohortName = 'NCG_high', 
  logscale = TRUE
)
dev.off()

### relative plot

get_TMB <- function(x) {
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  df <- x
  data_1 <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data_1)
}

cesc.csv  = get_TMB(x = CESC_maf)

cesc.csv$Barcode = substr(cesc.csv$Tumor_Sample_Barcode, 1, 12)
cesc.csv = cesc.csv %>% 
  filter(Barcode %in% life_data$Tumor_Sample_Barcode)
cesc.csv$Tumor_Sample_Barcode = cesc.csv$Barcode
cesc.csv = cesc.csv[-5]

newd = cesc.csv %>% 
  inner_join(life_data, by = "Tumor_Sample_Barcode")

newd$TMB = log2(newd$TMB + 1)

pdf(file="./Figure/step6_4.pdf", width = 8)
tiff(file = "./Figure/step6_4.tiff", 
     width = 10 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
newd %>% 
  ggplot(aes(x = Group, y = TMB)) +
  geom_boxplot(aes(fill = Group), width = 0.5, alpha = 0.4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        ) +
  plot.format +
  scale_fill_manual(values = c("HRG" = "#DDDDDD", "LRG" = "#E6B745")) +
  ggtitle("TMB") +
  labs(y = "Score")
dev.off()


all_data = all_data[-c(1:4)]
all_data = cbind(life_data, all_data)

full_data = cesc.csv %>% 
  as_tibble() %>% 
  inner_join(all_data, by="Tumor_Sample_Barcode")

subpddata = full_data[c(3,7)]
subgenedata <- full_data[, colnames(full_data) %in% cg]

new_da = cbind(subpddata, subgenedata)
save(new_da, file = "data/new_da.Rdata")
load("data/new_da.Rdata")

new_da$TMB = log(new_da$TMB + 1)
new_da[c(3:5)] = log(new_da[c(3:5)] + 1)

library(ggstatsplot)
pdf(file = "./Figure/step4_1.pdf", onefile = F)
ggscatterstats(data = new_da,
               y = TTC3,
               x = TMB,
               grouping.var = Group,
               #type = "pearson",
               results.subtitle = F,
               centrality.para = "mean",
               margins = "both",
               marginal.type = "densigram",
               title = "Relationship between TTC3 and TMB")

dev.off()

pdf(file = "./Figure/step4_2.pdf", onefile = F)
ggscatterstats(data = new_da,
               y = MGST1,
               x = TMB,
               grouping.var = Group,
               #type = "pearson",
               centrality.para = "mean",
               margins = "both",
               marginal.type = "densigram",
               title = "Relationship between MGST1 and TMB")
dev.off()

pdf(file = "./Figure/step4_3.pdf", onefile = F)
ggscatterstats(data = new_da,
               y = ACTR3,
               x = TMB,
               grouping.var = Group,
               #type = "pearson",
               centrality.para = "mean",
               margins = "both",
               marginal.type = "densigram",
               title = "Relationship between ACTR3 and TMB")
dev.off()

H = c("TTN", "PIK3CA", "KMT2C", "MUC16", "SYNE1", "DST",
      "EP300", "LRP1B", "NAV3", "ADGRV1", "FLG", "KMT2D",
      "DMD", "FBXW7", "RYR2", "USH2A", "DNAH2", "ERBB3",
      "NEB", "PTEN")

L = c("TTN", "PIK3CA", "KMT2C", "MUC16", "DMD", "FBXW7",
      "RYR2", "LRP2", "OBSCN", "EP300", "FAT1", "HUWE1",
      "USH2A", "ADGRV1", "DST", "MUC17", "FLG", "FLNA",
      "PCLO", "KMT2D")

H = data.frame(gene = H)
L = data.frame(gene = L)
cg = data.frame(gene = cg)

com = H %>% 
  inner_join(L, by = "gene")
L %>% 
  filter(!gene %in% com$gene) %>% 
  filter(gene %in% markgene2)

cg %>% 
  filter(gene %in% markgene2)
