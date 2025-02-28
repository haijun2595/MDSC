# 加载所需库
library(AUCell)
library(Seurat)
library(SingleCellExperiment)
suppressPackageStartupMessages({
  library("plyr")
  library("dplyr")
  library("ggpubr")
  library("ggsci")
  library("reshape2")
  library("survival")
  library("survminer")
  library("data.table")
})
library(GSVA)
library(ggbeeswarm)

input_dir <- "/input_TCGA/"
output_dir <- "/out"

cancer_types <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)

geneSets = list(
  "Fibroblast Score" = c("COL1A1", "ACTA2", "POSTN", "COL6A1"),
    'PMN_MDSC Score' = c("AQP9", "BCL2A1", "C5AR1", "CXCL1", "CXCR1", "CXCR2", "ADGRE2",
                       "FCGR3B", "FPR1", "FPR2", "ADGRG3", "HCK", "ICAM1", "IL1B", "IL1R2",
                       "CXCL8", "LILRB2", "LILRB3", "LYN", "MEFV", "NCF2", "OSM", "PLAUR",
                       "PTAFR", "S100A12", "S100A8", "S100A9", "SLC11A1", "SOD2", "TREM1", 
                       "TIMP2", "STAT3", "STAT6", "IRF1", "PTGS2", 
                       "IL4R", "MCEMP1", "NFKBIA", "TGFB1", "VEGFA"),
  'Exhaustion Score' = c("PDCD1", "LAYN", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1",
                         "CD3D", "CD3E","CD3G","CD8A","CD4","CD8B"),
'Antigen_Presentation Score'=c("ACTR1A", "ACTR1B", "AP1B1", "AP1M2", "AP2B1", "AP2M1", "CAPZB", "CD74", "CTSF", "CTSH", "CTSO", "DCTN2", "DNM3", "DYNC1I1", "DYNC1I2", "DYNC1LI1", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLADQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "KIF26A", "KIF2A", "KIF3B", "KIF3C", "KIF4B", "KIF5B", "KIFAP3", "KLC1", "KLC4", "SAR1B", "SEC13", "SPTBN2", "TUBA1A", "TUBA1B", "TUBA3C", "TUBA3D", "TUBA4A", "TUBA4B", "TUBA8", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B", "TUBB8", "TUBB8B"),
'Immunosuppressive Score'=c("IL1B", "IDO1", "CYBB", "NOS2", "ARG1", "TGFB1", "IL6", "IL10", "CD274", "VEGFA", "IL2RA", "CTLA4", "EBI3", "NT5E", "ENTPD1", "LILRB4", "LILRB2", "TDO2", "CD276", "VSIR", "LGALS9", "TNFRSF14", "PDCD1", "PDCD1LG2", "TREM1", "TREM2", "FOLR2"),
'EMT Score'=c("ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF", "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CCN1", "CCN2", "CD44", "CD59", "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COLGALT1", "COMP", "COPA", "CRLF1", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CXCL8", "DAB2", "DCN", "DKK1", "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3", "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1", "GPC1", "GPX7", "GREM1", "GREM1", "GREM1", "HTRA1", "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6", "INHBA", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15", "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", "MMP1", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9", "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "P3H1", "PCOLCE", "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1", "PLOD2", "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2", "PRSS2", "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1", "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1", "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3", "SNAI2", "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI", "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1", "VCAN", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A"
)
)


for (cancer_type in cancer_types) {
  cancer_rdata <- file.path(input_dir, cancer_type, paste0("merged_data_", cancer_type, ".rdata"))
  
  scRNA <- readRDS(cancer_rdata)
  
  scRNA <- scRNA[, -(1:3)]
  
  colnames(scRNA)[1] <- 'geneNames'
  
  sample_ids <- scRNA[, 1]  
  expr_matrix <- scRNA[, -1]  
  
  transposed_matrix <- t(expr_matrix)
  transposed_matrix <- as.data.frame(transposed_matrix)
  
  colnames(transposed_matrix) <- sample_ids
  
  transposed_matrix$cancer.type.abbreviation <- cancer_type 
  head(transposed_matrix)
  
  se.c <- transposed_matrix
  tmp <- gsva(as.matrix(se.c[, -c(ncol(se.c))]), geneSets, method="ssgsea", kcdf="Gaussian", ssgsea.norm=T, parallel.sz=40)
  
  bdf <- data.frame(t(tmp))
  bdf$sample <- se.c$sample_ids
  bdf$cancer.type.abbreviation <- rep(cancer_type, length(sample_ids))  
  

  colnames(bdf)[colnames(bdf) == "Fibroblast Score"] <- "Fibroblast.Score"
  colnames(bdf)[colnames(bdf) == "PMN_MDSC Score"] <- "PMN_MDSC.Score"
  colnames(bdf)[colnames(bdf) == "Exhaustion Score"] <- "Exhaustion.Score"
  colnames(bdf)[colnames(bdf) == "Antigen_Presentation Score"] <- "Antigen_Presentation.Score"
  colnames(bdf)[colnames(bdf) == "Immunosuppressive Score"] <- "Immunosuppressive.Score"
  colnames(bdf)[colnames(bdf) == "EMT Score"] <- "EMT.Score"
  
  p2 <- ggplot(bdf, aes(x=PMN_MDSC.Score, y=Fibroblast.Score)) +
    geom_point(shape=16, size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm', formula=y~x) +
    theme_classic2() +
    theme(legend.position="none") +
    ggtitle(cancer_type)
  ggsave(filename = file.path(output_dir, paste0(cancer_type, "_PMN_MDSC_Fibroblast.pdf")), p2, width = 4.5, height = 4)
  


  p4 <- ggplot(bdf, aes(x=PMN_MDSC.Score, y=Exhaustion.Score)) +
    geom_point(shape=16, size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm', formula=y~x) +
    theme_classic2() +
    theme(legend.position="none") +
    ggtitle(cancer_type)
  ggsave(filename = file.path(output_dir, paste0(cancer_type, "_PMN_MDSC_Exhaustion.pdf")), p4, width = 4.5, height = 4)
  
  p12 <- ggplot(bdf, aes(x=PMN_MDSC.Score, y=Antigen_Presentation.Score)) +
    geom_point(shape=16, size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm', formula=y~x) +
    theme_classic2() +
    theme(legend.position="none") +
    ggtitle(cancer_type)
  ggsave(filename = file.path(output_dir, paste0(cancer_type, "_PMN_MDSC_Antigen_Presentation.Score.pdf")), p12, width = 4.5, height = 4)
  
  p13 <- ggplot(bdf, aes(x=PMN_MDSC.Score, y=Immunosuppressive.Score)) +
    geom_point(shape=16, size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm', formula=y~x) +
    theme_classic2() +
    theme(legend.position="none") +
    ggtitle(cancer_type)
  ggsave(filename = file.path(output_dir, paste0(cancer_type, "_PMN_MDSC_Immunosuppressive.Score.pdf")), p13, width = 4.5, height = 4)
 
  p17 <- ggplot(bdf, aes(x=PMN_MDSC.Score, y=EMT.Score)) +
    geom_point(shape=16, size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm', formula=y~x) +
    theme_classic2() +
    theme(legend.position="none") +
    ggtitle(cancer_type)
  ggsave(filename = file.path(output_dir, paste0(cancer_type, "_PMN_MDSC_EMT.Score.pdf")), p17, width = 4.5, height = 4) 
 
}