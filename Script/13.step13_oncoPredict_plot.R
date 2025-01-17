rm(list = ls())

library(oncoPredict)

trainingExprData=readRDS(file='data/DataFiles/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds')
dim(trainingExprData) #
trainingExprData[1:4,1:4]

#IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already, and the calcPhenotype function Paul #has will do a power transformation (I assumed it would be better to not have both transformations)
trainingPtype = readRDS(file='data/DataFiles/DataFiles/Training Data/CTRP2_Res.rds')
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]

load("data/tpm_step5_all_survial_data.Rdata")

all_data[1:4,1:4]
data = all_data[-c(1:4)]
data = as.data.frame(t(data))

data = as.matrix(data)

save(trainingExprData, trainingPtype, data, file = "data/oncop.Rdata")

batchCorrect<-"eb"
powerTransformPhenotype<-TRUE
removeLowVaryingGenes<-0.2
removeLowVaringGenesFrom<-"homogenizeData"
minNumSamples=10
selection<- 1
printOutput=TRUE
pcr=FALSE
report_pc=FALSE
cc=FALSE
rsq=FALSE
percent=80

calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=data,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              pcr=pcr,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              report_pc=report_pc,
              percent=percent,
              rsq=rsq)


library(data.table)
testPtype <- fread('./data/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]

column_means <- colMeans(testPtype[-1], na.rm = TRUE)
top_20_min_columns <- names(sort(column_means)[1:20])
top_20_min_means <- column_means[top_20_min_columns]
top_20_min_means

testPtypemin20 = testPtype[, top_20_min_columns]
testPtypemin20 = log(testPtypemin20)
testPtypemin20$Barcode = testPtype$V1

load("data/step6_all_data.Rdata")
testPtypemin20$group = all_data$group
testPtypemin20$group = ifelse(testPtypemin20$group == "NCG_low", "LRG", "HRG")

rownames(testPtypemin20) = testPtypemin20$Barcode
testPtypemin20 = testPtypemin20[-21]
#testPtypemin20$group = factor(testPtypemin20$group, levels = c("NCG_low", "NCG_high"))

y_positions <- testPtypemin20 %>%
  pivot_longer(cols = -group, names_to = "chemic_type", values_to = "log(IC50_Value)") %>%
  group_by(chemic_type) %>%
  summarise(max_value = max(`log(IC50_Value)`) * 1.3) %>% 
  pull(max_value)

colnames(testPtypemin20) = c("Panobinostat(LBH589)", "Leptomycin B", "Dinaciclib",
                             "Ouabain", "CR-1-31B", "SB-743921", "Omacetaxine Mepesuccinate",
                             "BI-2536", "Doxorubicin:Navitoclax", "SR-II-138A",
                             "Flavopiridol(Alvocidib)", "Doxorubicin", "GSK461364",
                             "Narciclasine", "Bardoxolone Methyl", "AZD7762",
                             "Obatoclax(GX15-070)", "Brefeldin A(BFA)", "RSL3",
                             "Cucurbitacin I", "group")

color_labels <- function(x) {
  ifelse(x %in% c("Panobinostat(LBH589)", "Leptomycin B", "Ouabain",
                  "CR-1-31B", "SB-743921", "Omacetaxine Mepesuccinate",
                  "SR-II-138A", "Flavopiridol(Alvocidib)", "Doxorubicin",
                  "GSK461364", "Narciclasine", "Bardoxolone Methyl",
                  "AZD7762", "Obatoclax(GX15-070)"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

pdf(file="./Figure/step5_4.pdf", width = 15)
tiff(file = "./Figure/step5_4.tiff", 
     width = 12 * 300,  
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
testPtypemin20 %>%
  pivot_longer(cols = -group, names_to = "chemic_type", values_to = "log(IC50_Value)") %>%
  ggplot(aes(x = chemic_type, y = `log(IC50_Value)`)) +
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
  ggtitle("oncoPredict") +
  labs(y = "log(IC50)")
dev.off()

