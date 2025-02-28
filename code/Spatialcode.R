library(Seurat)
library(data.table)
library(tidyverse)
library(MCPcounter)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

#devtools::install_github("MarioniLab/DropletUtils")
library(DropletUtils)
library(tibble)
library(jsonlite)
library(stringr)
#install_github("ebecht/MCPcounter",ref="master", subdir="Source")

dir = "D:/New/"
ST_object <- Load10X_Spatial(data.dir = dir,                            
                             filename = "filtered_feature_bc_matrix.h5")

for (i in colnames((ST_object@images$slice1@coordinates))) {
  ST_object@images$slice1@coordinates[[i]] <- as.integer(ST_object@images$slice1@coordinates[[i]])
} 
SpatialDimPlot(ST_object,alpha = 0)

ST_object <- SCTransform(ST_object, assay = "Spatial", verbose = FALSE)
ST_object <- RunPCA(ST_object, assay = "SCT", verbose = FALSE) 
ST_object <- FindNeighbors(ST_object, reduction = "pca", dims = 1:30)
ST_object <- FindClusters(ST_object, verbose = FALSE) 
ST_object <- RunUMAP(ST_object, reduction = "pca", dims = 1:30 )
head(ST_object)
save(ST_object,file="ST1.Rdata")
#SpatialDimPlot(ST_object,alpha = 1)


markerList <- list(
  "Fibroblast" = c("COL1A1", "ACTA2", "POSTN", "COL6A1"),
  "Myeloid" = c("CD68","CD14", "C1QA", "LYZ","APOE"),
  "Epithelial_Tumor_cell"=c("EPCAM","KRT8","CDH16","CA9","ANGPTL4"
                 #"KRT18","ERBB2","NAPSA","KRT7",
                 ),
  "PMN_MDSC"= c("AQP9", "BCL2A1", "C5AR1", "CXCL1", "CXCR1", "CXCR2", "ADGRE2",
                "FCGR3B", "FPR1", "FPR2", "ADGRG3", "HCK", "ICAM1", "IL1B", "IL1R2",
                "CXCL8", "LILRB2", "LILRB3", "LYN", "MEFV", "NCF2", "OSM", "PLAUR",
                "PTAFR", "S100A12", "S100A8", "S100A9", "SLC11A1", "SOD2", "TREM1", 
                "TIMP2",  "STAT3", "STAT6", "IRF1", "PTGS2", 
                "IL4R", "MCEMP1", "NFKBIA", "TGFB1", "VEGFA"),
  "Tex"=c("TIGIT","PDCD1","CTLA4","HAVCR2","LAG3","LAYN","TOX","VSIR","BTLA","ENTPD1",
          "CD3D", "CD3E","CD3G","CD8A","CD4","CD8B"),
  "Immunosuppressive"=c("IL1B", "IDO1", "CYBB", "NOS2", "ARG1", "TGFB1", "IL6", "IL10", 
                        "CD274", "VEGFA", "IL2RA", "CTLA4", "EBI3", "NT5E", "ENTPD1", 
                        "LILRB4", "LILRB2", "TDO2", "CD276", "VSIR", "LGALS9", "TNFRSF14", 
                        "PDCD1", "PDCD1LG2", "TREM1", "TREM2", "FOLR2"),
  "EMT"=c("ABI3BP", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF", 
          "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CCN1", "CCN2", 
          "CD44", "CD59", "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", 
          "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", 
          "COL5A3", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COLGALT1", "COMP", 
          "COPA", "CRLF1", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CXCL8", "DAB2", 
          "DCN", "DKK1", "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", 
          "EMP3", "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", 
          "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", 
          "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1", 
          "GPC1", "GPX7", "GREM1", "GREM1", "GREM1", "HTRA1", "ID2", "IGFBP2", 
          "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6", "INHBA", "ITGA2", "ITGA5", 
          "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", 
          "LAMC1", "LAMC2", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15", 
          "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", 
          "MMP1", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9", 
          "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "P3H1", "PCOLCE", 
          "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1", "PLOD2", 
          "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2", "PRSS2", 
          "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1", "SCG2", "SDC1", 
          "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1", "SFRP4", "SGCB", 
          "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3", "SNAI2", "SNTB1", "SPARC", 
          "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI", "TGFBR3", "TGM2", 
          "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", 
          "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1", "VEGFA","VCAN", "VEGFC", 
          "VIM", "WIPF1", "WNT5A")
)


markerList2 <- markerList[c("Fibroblast","Epithelial_Tumor_cell",
                            "Myeloid", "PMN_MDSC","Tex","Immunosuppressive","EMT")]
for (i in names(markerList2)) {
  intersect_sign<-intersect(markerList2[[i]],rownames(ST_object))
  cell_types <- make.names(paste(i ,"GeneMean",sep = "_"))#
  cell_types_Seurat <- make.names(paste(i ,"Seurat",sep = "_"))  

  ST_object <- AddMetaData(ST_object,apply(as.matrix(ST_object[["SCT"]]@data[intersect_sign,]),2,mean),
                           col.name = cell_types  ) 
  # AddModuleScore #
  ST_object<-AddModuleScore(object=ST_object,features=list(intersect_sign),name=cell_types_Seurat)
}
head(ST_object)
view(ST_object@meta.data)
##########################################################
##########################################################
##########################################################


output_dir <- dir

for (cell_type in names(markerList2)) {
  feature_name <- make.names(paste(cell_type, "Seurat1", sep = "_"))
  
  if (feature_name %in% colnames(ST_object@meta.data)) {
    p <- SpatialPlot(ST_object,
                     features = feature_name,
                     image.alpha = 0,  
                     pt.size.factor = 1.7,
                     alpha = 1) +
      theme(
        text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
    
    output_file <- paste0(output_dir, cell_type, "_SpatialPlot.pdf")
    ggsave(output_file, plot = p, width = 8, height = 7)
    
    print(paste("Saved SpatialPlot for:", cell_type, "to", output_file))
  } else {
    print(paste("Warning: Feature", feature_name, "not found in ST_object@meta.data"))
  }
}
