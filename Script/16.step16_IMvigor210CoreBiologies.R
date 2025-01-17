rm(list = ls())

load("./Rdata/IMvigor210CoreBiologies.Rdata")

annoData = na.omit(annoData)

expr = exprSet %>% 
  rownames_to_column("entrez_id") %>% 
  as_tibble() %>% 
  inner_join(annoData, by = "entrez_id") %>% 
  as.data.frame()

expr = expr[-c(1, 351:354)]

nor_data = expr
nor_data$mean_expression <- rowMeans(nor_data[1:348])
nor_data <- nor_data[order(nor_data$symbol, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$symbol), ]
rownames(nor_data) = nor_data$symbol
nor_data = nor_data[1:348]

nor_data = nor_data[-1, ]
nor_data = as.data.frame(t(nor_data))

nor_data = nor_data %>% 
  rownames_to_column("sample")

all_data = phenoData %>% 
  rownames_to_column("sample") %>% 
  inner_join(nor_data, by = "sample")

rownames(all_data) = all_data$sample

table(all_data$Tissue)

all_data = filter(all_data, Tissue == "bladder")

all_data$GADH11111 = (3*(all_data$PERP + all_data$ANXA1)) / (2*(all_data$ARGLU1 + all_data$S100A6 + all_data$`HLA-E`))
cg = c("GADH11111")

data = all_data[-c(1:26)]
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

life_data = all_data[1:26]
life_data$group = as.character(group.gsva)
life_data$group = factor(life_data$group, levels = c('LRG','HRG'))
life_data$os = life_data$os * 12

cox_fit <- coxph(Surv(os, censOS) ~ group, data = life_data)
hr <- round(summary(cox_fit)$conf.int[1],3)
hr_lower <- round(summary(cox_fit)$conf.int[3], 3)
hr_upper <- round(summary(cox_fit)$conf.int[4], 3)

pdf(file = "./Figure/step7_7.pdf", onefile = F, width = 8)
tiff(file = "./Figure/step7_7.tiff", 
     width = 10 * 300,  
     height = 8 * 300,  
     res = 300, 
     compression = "lzw")
survfit2(Surv(os, censOS) ~ group, data = life_data) %>% 
  ggsurvfit(linewidth = 0.8) +
  add_pvalue(location  = "annotation", x = 270, y = 0.88 , hjust = 1,
             size = 4, caption = "Stratified log-rank {p.value}") +
  annotate("text",    
           x = 270,
           y = 0.95,
           label = paste("Hazard ratio,", round(hr,2), "(95% CI,", round(hr_lower,2), "-", round(hr_upper,2), ")", sep = ""),
           size = 4,
           hjust = 1
  ) +
  annotate("text",   
           x = 200,
           y = 0.55,
           label = "LRG",
           size = 4.5,
           hjust = 0,
  ) +  
  annotate("text",   
           x = 100,
           y = 0.4,
           label = "HRG",
           size = 4.5,
           hjust = 1  
  ) +  
  labs(title = "IMvigor210CoreBiologies-bladder(n=195)",
       x = "Months",
       y = "overall survival",
  ) +
  scale_x_continuous(breaks = seq(0, 290, 30), expand = c(0.03, 0)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  scale_color_manual(values = c('#2D5662', '#F7931D')) +
  scale_fill_manual(values = c('#2D5662', '#F7931D')) +
  theme_classic() +
  theme(axis.text = element_text(size = 14,color = "black"), 
        axis.title.y = element_text(size = 14,color = "black", vjust = 1),
        axis.title.x = element_text(size = 14,color = "black"),
        panel.grid.major.y =  element_line(color = "#DCDDDF"),
        legend.position = "none",
        plot.title=element_text(color="black",size=30, face = "bold")
        )
dev.off()

###
all_data$group = as.character(group.gsva)

Cellratio = prop.table(table(all_data$binaryResponse, all_data$group), margin = 2)
Cellratio = as.data.frame(Cellratio)

colnames(Cellratio) = c("Class", "Group", "Freq")

pdf(file = "./Figure/step7_8.pdf", width = 4)
ggplot(Cellratio, aes(Group, Freq, fill = Class)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Response", x = "", y = "Estimated Proportion") + 
  theme_bw() +
  theme(legend.position = 'right') +
  plot.format + 
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = c("CR/PR" = "#DDDDDD", "SD/PD" = "#E6B745"))
dev.off()
