
library(GSVA)     
library(ggplot2)  
library(reshape2) 
library(ggpubr)  

immune <- read.table("CIBERSORT-Results.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

immune <- immune[immune$`P-value` < 0.05, ]

valid_samples <- rownames(immune)

expr <- as.matrix(read.table("uniq.symbol.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

expr <- expr[, colnames(expr) %in% valid_samples]

dim(expr)


gene_sets_file <- "geneset.txt" 
gene_sets <- read.table(gene_sets_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


gene_set_list <- lapply(gene_sets, function(col) {
  genes <- na.omit(col)  
  return(genes)
})

names(gene_set_list) <- colnames(gene_sets)

head(gene_set_list)
set.seed(123456)
gsva_scores <- gsva(expr, gene_set_list, method = "ssgsea", verbose = FALSE)

head(gsva_scores)

mdsc_score <- gsva_scores["PMN_MDSC", ]
median_mdsc_score <- median(mdsc_score)
mdsc_group <- ifelse(mdsc_score > median_mdsc_score, "high", "low")

mdsc_group_info <- data.frame(Sample = colnames(expr),
                              PMN_MDSC_Score = mdsc_score, 
                              PMN_MDSC_Level = mdsc_group)

merged_data_all <- data.frame(Sample = colnames(expr),
                              PMN_MDSC_Level = mdsc_group_info$PMN_MDSC_Level)

merged_data_all$PMN_MDSC_Level <- factor(merged_data_all$PMN_MDSC_Level, levels = c("low","high"))

for (gene_set_name in names(gene_set_list)) {
  merged_data_all[gene_set_name] <- gsva_scores[gene_set_name, ]
}


gene_set_name <- names(gene_set_list)[1]
low_group <- merged_data_all[merged_data_all$PMN_MDSC_Level == "low", gene_set_name]
high_group <- merged_data_all[merged_data_all$PMN_MDSC_Level == "high", gene_set_name]

shapiro_low <- shapiro.test(low_group)
shapiro_high <- shapiro.test(high_group)

print(shapiro_low)
print(shapiro_high)
library(car)
levene_test <- leveneTest(as.formula(paste(gene_set_name, "~ PMN_MDSC_Level")), data = merged_data_all)
print(levene_test)




library(ggplot2)
library(rlang) 

library(ggplot2)
library(ggpubr)
for (gene_set_name in names(gene_set_list)) {
  ridge_plot <- ggplot(merged_data_all, aes(x = !!sym(gene_set_name), fill = PMN_MDSC_Level)) + 
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("high" = "#6e336e", "low" = "#008a5d")) +
    labs(title = paste(gene_set_name, "Score"),
         x = paste(gene_set_name, "Score")) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "top",
      panel.grid = element_blank(),  
      axis.title.y = element_blank() 
    )
  
  ggsave(paste0("Density_Plot_", gene_set_name, ".png"), 
         plot = ridge_plot, 
         width = 5, 
         height = 4)  
  

  boxplot_plot <- ggplot(merged_data_all, aes(x = PMN_MDSC_Level, y = !!sym(gene_set_name), fill = PMN_MDSC_Level)) +
    geom_boxplot(alpha = 0.6, width = 0.4) +
    stat_compare_means(aes(group = PMN_MDSC_Level), label = "p.signif", method = "wilcox.test",comparisons = list(c("low", "high"))) +
    scale_fill_manual(values = c("high" = "#6e336e", "low" = "#008a5d")) +
    labs(title = paste(gene_set_name, "Score"),
         x = "PMN_MDSC Level",
         y = paste(gene_set_name, "Score")) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "top",
      panel.grid = element_blank(),  
      axis.title.y = element_text(color = "black", size = 12), 
      axis.title.x = element_text(color = "black", size = 12),  
      axis.line = element_line(colour = "black"), 
      axis.ticks = element_line(colour = "black") 
    )
  

  ggsave(paste0("Boxplot_", gene_set_name, ".pdf"), 
         plot = boxplot_plot, 
         width = 4,  
         height = 4)  
}


