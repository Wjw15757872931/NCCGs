## IPS

rm(list = ls())

library(tidyverse)

pd = read_tsv("./data/TCIA-ClinicalData.tsv")

pd = pd[c(1, 51:54)]

load("./data/step5_all_survial_data.Rdata")

cg = c("TTC3", "MGST1", 'ACTR3')

data = all_data[-c(1:4)]
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
group.gsva <- ifelse( dat.gsva < 0, "LRG", "HRG")
table(group.gsva)

life_data = all_data[1:4]
life_data$group = as.character(group.gsva)

data = all_data[-c(1:4)]

all_data = cbind(life_data, data)
save(all_data, file = "data/step6_all_data.Rdata")

fd = all_data[4:5]
colnames(fd) = c("barcode", "group")

newd = fd %>% 
  inner_join(pd, by="barcode")

rownames(newd) = newd$barcode
newd = newd[-1]
#newd$group = factor(newd$group, levels = c("NCG_low", "NCG_high"))

y_positions <- newd %>%
  pivot_longer(cols = -group, names_to = "IPS", values_to = "Value") %>%
  group_by(IPS) %>%
  summarise(max_value = max(Value) * 1.3) %>% 
  pull(max_value)

colnames(newd) = c("group", "CTLA4(-)/PD-1(-)", "CTLA4(-)/PD-1(+)", "CTLA4(+)/PD-1(-)", "CTLA4(+)/PD-1(+)")

color_labels <- function(x) {
  ifelse(x %in% c("CTLA4(-)/PD-1(+)", "CTLA4(+)/PD-1(+)"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}


pdf(file="./Figure/step5_1.pdf", width = 10)
tiff(file = "./Figure/step5_1.tiff", 
     width = 8 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
newd %>% 
  pivot_longer(cols = -group, names_to = "IPS", values_to = "Value") %>%
  ggplot(aes(x = IPS, y = Value)) +
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
  ggtitle("TCIA-IPS") +
  labs(y = "Score")
dev.off()

### TIDE
rm(list = ls())
load("data/step6_all_data.Rdata")

data = all_data[-c(1:5)]
data = as.data.frame(t(data))

normalize <- t(apply(data, 1, function(x)x-(mean(x))))
write.table(normalize, file = "data/exp_tide.txt", sep = "\t", 
            row.names = T, col.names = TRUE, quote = FALSE)

tide = read.csv("data/tide_result.csv")
pd = all_data[c(1,5)]
colnames(pd) = c("Patient", "group")
tide = tide %>% 
  inner_join(pd, by = "Patient")

Cellratio = prop.table(table(tide$Responder, tide$group), margin = 2)
Cellratio = as.data.frame(Cellratio)

colnames(Cellratio) = c("Class", "Group", "Freq")
Cellratio$Group = ifelse(Cellratio$Group == "NCG_low", "LRG", "HRG")
#Cellratio$Group = factor(Cellratio$Group, levels = c("NCG_low", "NCG_high"))

pdf(file = "./Figure/step5_2.pdf", width = 4)
ggplot(Cellratio, aes(Group, Freq, fill = Class)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Responder", x = "", y = "Estimated Proportion") + 
  theme_bw() +
  theme(legend.position = 'right') +
  plot.format + 
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = c("False" = "#DDDDDD", "True" = "#E6B745"))
dev.off()

Cellratio = prop.table(table(tide$CTL.flag, tide$group), margin = 2)
Cellratio = as.data.frame(Cellratio)

colnames(Cellratio) = c("Class", "Group", "Freq")
Cellratio$Group = ifelse(Cellratio$Group == "NCG_low", "LRG", "HRG")
#Cellratio$Group = factor(Cellratio$Group, levels = c("NCG_low", "NCG_high"))

pdf(file = "./Figure/step5_5.pdf", width = 4)
ggplot(Cellratio, aes(Group, Freq, fill = Class)) + 
  geom_bar(stat = "identity") +
  labs(fill = "CTL.flag", x = "", y = "Estimated Proportion") + 
  theme_bw() +
  theme(legend.position = 'right') +
  plot.format + 
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = c("False" = "#DDDDDD", "True" = "#E6B745"))
dev.off()


tide_1 = tide[c(4, 11:12, 16)]
tide_1$group = ifelse(tide_1$group == "NCG_low", "LRG", "HRG")

y_positions <- tide_1 %>%
  pivot_longer(cols = -group, names_to = "type", values_to = "Score") %>%
  group_by(type) %>%
  summarise(max_value = max(Score) * 1.3) %>% 
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("TIDE", "Dysfunction", "Exclusion"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

#tide_1$group = factor(tide_1$group, levels = c("NCG_low", "NCG_high"))
pdf(file="./Figure/step5_3.pdf", width = 8)
tiff(file = "./Figure/step5_3.tiff", 
     width = 8 * 300,   
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")

tide_1 %>% 
  pivot_longer(cols = -group, names_to = "type", values_to = "Score") %>%
  ggplot(aes(x = type, y = Score)) +
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
  ggtitle("TIDE") +
  labs(y = "Score")
dev.off()
