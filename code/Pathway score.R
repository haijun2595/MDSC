library(Seurat)
library(tidyverse)
library(pheatmap)
#gene_signatures_path <- "gene signature.txt"
load("/Neu_sub_Dele_harmony.Rdata")
gene_signatures_path <- "gene signature.txt"


gene_data <- read.table(gene_signatures_path, header = TRUE, sep = "\t", check.names = FALSE)
library(Seurat)
library(tidyverse)
library(pheatmap)

marker_list <- list()

for (i in 1:ncol(gene_data)) {
  gene_set_name <- colnames(gene_data)[i]
  genes <- na.omit(gene_data[[i]])
  marker_list[[gene_set_name]] <- genes
}


print(marker_list)
Idents(Neu_sub_Dele_harmony) <- Neu_sub_Dele_harmony$subcelltype
Neu_sub_Dele_harmony <- AddModuleScore(Neu_sub_Dele_harmony, 
                         features = marker_list, 
                         ctrl = 5, 
                         name = "FunctionScore")
for (i in 1:length(marker_list)) {
  colnames(Neu_sub_Dele_harmony@meta.data)[colnames(Neu_sub_Dele_harmony@meta.data) == paste0("FunctionScore", i)] <- names(marker_list)[i]
}


colnames(Neu_sub_Dele_harmony@meta.data)
head(Neu_sub_Dele_harmony@meta.data)
library(Seurat)
library(tidyverse)
library(pheatmap)


Idents(Neu_sub_Dele_harmony) <- Neu_sub_Dele_harmony$subcelltype

Function <- c("Migration", "Proliferation", "IFN response","Stemness","Chemotaxis"
              #"Neutrophil Maturation","Neutrophil Aging"
              )
ROS <- c("ROS generation", "ROS quenching")

MDSC_Function <- c("MDSC function")
NADPH <- c("NADPH oxidation","NADPH oxidase complex"
           #"NADPH oxidase H2O2 forming activity"
           )
Metabolism <- c("OXPHOS", "Glycolysis", "FAO")
Immunosuppressive <- c("Immunosuppressive"
                       #"MDSCImmunosuppression"
                       )
VEGF <- c("VEGF")

MarkerNameVector <- c(Function, ROS,MDSC_Function,NADPH,
                      Metabolism, Immunosuppressive,VEGF)

View(Neu_sub_Dele_harmony@meta.data)

unique_clusters <- sort(unique(Neu_sub_Dele_harmony$subcelltype))
table(Neu_sub_Dele_harmony$subcelltype)
FunctionScoreMatrix <- matrix(0, 
                              ncol = length(unique_clusters), 
                              nrow = length(MarkerNameVector)) 
colnames(FunctionScoreMatrix) <- unique_clusters
rownames(FunctionScoreMatrix) <- MarkerNameVector


for (ci in 1:ncol(FunctionScoreMatrix)) {
  for (ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(Neu_sub_Dele_harmony@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[Neu_sub_Dele_harmony$subcelltype == unique_clusters[ci]], na.rm = TRUE)
    FunctionScoreMatrix[ri, ci] <- fv
  }
}


FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scale))




my.breaks <- seq(-1, 1, by = 0.1)
my.colors <- colorRampPalette(colors = c("#6DCCFD", "white", "#FD9AA0"))(length(my.breaks))
#my.colors <- colorRampPalette(colors = c("#6DCCFD",  "#EEEEEE", "#FFAAAA", "#E15759"))(length(my.breaks))
#my.colors <- colorRampPalette(colors = c("#6DCCFD", "white", "#FD9AA0"))(length(my.breaks))
#my.colors <- colorRampPalette(colors = c("#6DCCFD",  "#EEEEEE", "#FFAAAA", "#E15759"))(length(my.breaks))
my.colors <- colorRampPalette(colors = c("#008bd0", "#A8D5E9","#eeeeee", "#FFBB78", "#ffa61d"))(length(my.breaks))
my.colors <- colorRampPalette(colors = c("#C2D95E","#57ab81", "white", "#ff9600"))(length(my.breaks))

signatureType_row <- data.frame(Signature.type = c(
  rep("Function", length(Function)),
  rep("ROS", length(ROS)),
  #rep("Chemotaxis", length(Chemotaxis)),
  rep("MDSC_Function", length(MDSC_Function)),
  rep("NADPH", length(NADPH)),
  rep("Metabolism", length(Metabolism)),
  #rep("Checkpoint", length(Checkpoint)),
  rep("Immunosuppressive", length(Immunosuppressive)),
  #rep("iNOS", length(iNOS)),
  #rep("Exhaustion", length(Exhaustion)),
  rep("VEGF", length(VEGF))
  ))

rownames(signatureType_row) <- MarkerNameVector

unique_clusters <- levels(Neu_sub_Dele_harmony$subcelltype)

colnames(FunctionScoreMatrix) <- unique_clusters


pdf(file="heatmap.pdf",width=12,height=4)
pheatmap(t(FunctionScoreMatrix),
         show_colnames = TRUE,
         show_rownames = TRUE,
         #annotation_row = signatureType_row,
         annotation_col = signatureType_row,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         #gaps_row = c(5, 6,9,10, 12),  
         gaps_col = c(5,7,8,10,13,14),   
         breaks = my.breaks,
         color = my.colors,
         border_color = "black",  
         fontsize = 8,
         fontsize_row = 10,        
         fontsize_col = 10,      
         angle_col = "45",       
         width = 8,              
         height = 6,           
         annotation_colors = list(
           Signature.type = c("Function" = "#3969AC", 
                              "ROS" = "lightgreen", 
                              "MDSC_Function"="#BDBADB",
                              "NADPH"="#D65190",
                              "Metabolism" = "#F28E2B", 
                              "Immunosuppressive" = "#AA40FC",
                              "VEGF"="#C2D95E"))
)

dev.off()




pdf(file="heatmap2.pdf",width=6,height=9)
pheatmap(FunctionScoreMatrix,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_row = signatureType_row,
         #annotation_col = signatureType_row,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row = c(7,9,10,13,16,18), 
         #gaps_col = c(7,9,10,13,16,18),
         breaks = my.breaks,
         color = my.colors,
         border_color = "black",  
         fontsize = 8,
         fontsize_row = 10,     
         fontsize_col = 10,      
         angle_col = "45",        
         width = 8,             
         height = 6,             
         annotation_colors = list(
           Signature.type = c("Function" = "#3969AC", 
                              "ROS" = "lightgreen", 
                              "MDSC_Function"="#BDBADB",
                              "NADPH"="#D65190",
                              "Metabolism" = "#F28E2B", 
                              "Immunosuppressive" = "#AA40FC",
                              "VEGF"="#C2D95E"))
)
dev.off()
