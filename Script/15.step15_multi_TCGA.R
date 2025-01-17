### multi TCGA HOT PLOT

rm(list = ls())

load("data/sn_sc_degs.Rdata")
scf1 = diff_DEG_sc %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25)

sc_sn_genes = diff_DEG_sn %>% 
  rownames_to_column("symbol") %>% 
  as_tibble() %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(pct.1 > 0.25) %>% 
  filter(pct.2 > 0.25) %>% 
  filter(symbol %in% scf1$rowname) %>% 
  dplyr::select(symbol) %>% 
  as.data.frame()

sc_sn_genes$'TCGA-CESC' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", "P4HB", "TTC3",
                                                         "TPM3", "TPT1", "TM9SF2",
                                                         "ACTR3", "MGST1") ~ "1",
                                  TRUE ~ "0"
                               )
sc_sn_genes$'HPV+ TCGA-CESC' = case_when(sc_sn_genes$symbol %in% c("TTC3", "DNAJC3", "ACTR3", "MGST1") ~ "1",
                                         sc_sn_genes$symbol == "TMSB4X" ~ "-1",
                                         TRUE ~ "0"
)

sc_sn_genes$'TCGA-HNSC' = case_when(sc_sn_genes$symbol %in% c("TPT1", "TPM3", "TMBIM6",
                                                              "TM9SF2", "PERP", "PDIA3",
                                                              "PABPC1", "P4HB", "GAPDH",
                                                              "CANX", "CALR", "ARPC2",
                                                              "ACTR3" ) ~ "1",
                                         TRUE ~ "0"
)

sc_sn_genes$'TCGA-BLCA' = case_when(sc_sn_genes$symbol %in% c("PERP", "ANXA1") ~ "1",
                                    sc_sn_genes$symbol %in% c("S100A6", "HLA-E", "ARGLU1") ~ "-1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-BRCA' = case_when(
                                    sc_sn_genes$symbol %in% c("CD74", "B4GALT1", "ARGLU1") ~ "-1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-COAD' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", "PDIA3", "HLA-E") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-ESCA' = case_when(sc_sn_genes$symbol %in% c("TPT1", "PPP1CB") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-KIRC' = case_when(sc_sn_genes$symbol %in% c("ARGLU1") ~ "1",
                                    sc_sn_genes$symbol %in% c("TSC22D1", "TMBIM6", "TM9SF2",
                                                              "RBM47", "PPP1CB", "HSD17B4",
                                                              "HLA-E", "ERGIC1", "CANX", "ATP1B1") ~ "-1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-KIRP' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", 'PABPC1', "IFI16",
                                                              "FBLN1", "ACTR3") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-LIHC' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", "VGLL4", "TPM3",
                                                              "TMSB4X", "TMSB10", "RBM47",
                                                              "PPP1CB", "PERP", "PABPC1",
                                                              "LGALS3", "GFPT1", "GAPDH",
                                                              "ARPC2", "ACTR3") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-LUAD' = case_when(sc_sn_genes$symbol %in% c("GAPDH") ~ "1",
                                         sc_sn_genes$symbol == "DHRS3" ~ "-1",
                                         TRUE ~ "0"
)

sc_sn_genes$'TCGA-LUSC' = case_when(sc_sn_genes$symbol %in% c("TSC22D1", "TM9SF2", "RNF213",
                                                              "HSD17B4", "HLA-E", "ERGIC1",
                                                              "DNAJC3", "DHRS3", "CD74",
                                                              "CANX", "B4GALT1", "ARPC2") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-OV' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", "TPT1", "ELF3") ~ "1",
                                    sc_sn_genes$symbol %in% c("B4GALT1") ~ "-1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-PAAD' = case_when(sc_sn_genes$symbol %in% c("ZFP36L1", "VGLL4", "TMSB10",
                                                              "PPP1CB", "PERP", "PDIA3",
                                                              "PABPC1", "LGALS3", "GAPDH",
                                                              "ELF3", "CANX", "B4GALT1") ~ "1",
                                    TRUE ~ "0"
)

sc_sn_genes$'TCGA-STAD' = case_when(sc_sn_genes$symbol %in% c("HSPB1", "FBLN1") ~ "1",
                                  sc_sn_genes$symbol %in% c("ARGLU1") ~ "-1",
                                  TRUE ~ "0"
)

sc_sn_genes$'TCGA-UCEC' = case_when(sc_sn_genes$symbol %in% c("RNF213", "PABPC1", "MGST1", "ARPC2", "ACTR3") ~ "1",
                                    sc_sn_genes$symbol %in% c("B4GALT1") ~ "-1",
                                    TRUE ~ "0"
)

rownames(sc_sn_genes) = sc_sn_genes$symbol
sc_sn_genes = sc_sn_genes[-1]

library(pheatmap)

sc_sn_genes[] <- lapply(sc_sn_genes, as.numeric)
sc_sn_genes = as.data.frame(t(sc_sn_genes))

# Replace numeric values with text
#text_matrix <- ifelse(sc_sn_genes == -1, "HR < 1",
#                      ifelse(sc_sn_genes == 0, "No effect", "HR > 1"))

pdf(file="./Figure/step7_6.pdf", width = 30, height = 12)
pheatmap(
  sc_sn_genes,
  cluster_rows = F,
  cluster_cols = T,
  color = c("#3674a2", "#FFFFFF", "#df5734"), 
  legend_breaks = c(-1, 0, 1),       
  legend_labels = c("HR < 1", "No effect", "HR > 1"), 
  #display_numbers = text_matrix,     
  #number_color = "black",
  cellwidth = 40, cellheight = 40,
  fontsize = 15
)
dev.off()
