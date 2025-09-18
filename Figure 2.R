rm(list=ls())
setwd(getwd())
library(pheatmap)
library(tidyverse)
library(ComplexHeatmap)
library(magick)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 2B
load('./Data/Figure 2B.rda')

png("./Results/Figure 2B.Pattern1a.1.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_1a_1
dev.off()

png("./Results/Figure 2B.Pattern1a.2.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_1a_2
dev.off()

png("./Results/Figure 2B.Pattern1b.1.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_1b_1
dev.off()

png("./Results/Figure 2B.Pattern1b.2.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_1b_2
dev.off()

png("./Results/Figure 2B.Pattern2.1.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_2_1
dev.off()

png("./Results/Figure 2B.Pattern2.2.png", width = 4, height = 8, res = 300, units = 'in')
Pattern_2_2
dev.off()


# Figure 2C
# P20
load('./Data/Figure 2C.P20.rda')

png("./Results/Figure 2C.P20_PhyloTree.png", width = 4, height = 8, res = 300, units = 'in')
PhyloTree
dev.off()

pheatmap(plot_mat, scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), breaks = seq(0.9,1.1,by=0.004),angle_col = "90", 
         cluster_cols = F, gaps_row = NULL, gaps_col = chromosome_info, annotation_colors = annotation_colors, fontsize = 8,
         annotation_row = annotation_row,
         annotation_col = annotation_col, 
         cluster_rows = F, show_colnames = F, show_rownames = F, border_color = NA,
         width = 8,
         height = 4,
         cutree_col = 2, 
         file = paste0("./Results/Figure 2C.P20_CNA_heatmap.png"))


# P3
load('./Data/Figure 2C.P3.rda')

png("./Results/Figure 2C.P3_PhyloTree.png", width = 4, height = 8, res = 300, units = 'in')
PhyloTree
dev.off()

pheatmap(plot_mat, scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), breaks = seq(0.9,1.1,by=0.004),angle_col = "90", 
         cluster_cols = F, gaps_row = NULL, gaps_col = chromosome_info, annotation_colors = annotation_colors, fontsize = 8,
         annotation_row = annotation_row,
         annotation_col = annotation_col, 
         cluster_rows = F, show_colnames = F, show_rownames = F, border_color = NA,
         width = 8,
         height = 4,
         cutree_col = 2, 
         file = paste0("./Results/Figure 2C.P3_CNA_heatmap.png"))


# P23
load('./Data/Figure 2C.P23.rda')

png("./Results/Figure 2C.P23_PhyloTree.png", width = 4, height = 8, res = 300, units = 'in')
PhyloTree
dev.off()

pheatmap(plot_mat, scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), breaks = seq(0.9,1.1,by=0.004),angle_col = "90", 
         cluster_cols = F, gaps_row = NULL, gaps_col = chromosome_info, annotation_colors = annotation_colors, fontsize = 8,
         annotation_row = annotation_row,
         annotation_col = annotation_col, 
         cluster_rows = F, show_colnames = F, show_rownames = F, border_color = NA,
         width = 8,
         height = 4,
         cutree_col = 2, 
         file = paste0("./Results/Figure 2C.P23_CNA_heatmap.png"))

# Figure 2C spatial plots were created by using 'Figure 2C.py'


# Figure 2D
load('./Data/Figure 2D.rda')

col = c('Missense_Mutation' = '#377eb8', 'Truncating' = '#ff7f00', 'In_Frame_InDel' = '#f781bf')
p1 <- oncoPrint(mat1, alter_fun = alter_fun, 
                col = col, 
                column_title = 'Pattern 1a', 
                top_annotation = NULL,
                heatmap_legend_param = list(title = 'Alternations', at = c('Missense_Mutation', 'Truncating', 'In_Frame_InDel'), labels = c('Missense', 'Truncating', 'In_Frame_InDel')),
                show_column_names = TRUE,
                show_row_names = T,
                pct_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12)) 

p2 <- oncoPrint(mat2, alter_fun = alter_fun, 
                col = col, 
                column_title = 'Pattern 1b',
                top_annotation = NULL,
                heatmap_legend_param = list(title = 'Alternations', at = c('Missense_Mutation', 'Truncating', 'In_Frame_InDel'), labels = c('Missense', 'Truncating', 'In_Frame_InDel')),
                show_column_names = TRUE,
                show_row_names = T,
                pct_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12)) 

p3 <- oncoPrint(mat3, alter_fun = alter_fun, 
                col = col, 
                column_title = 'Pattern 2',
                top_annotation = NULL,
                heatmap_legend_param = list(title = 'Alternations', at = c('Missense_Mutation', 'Truncating', 'In_Frame_InDel'), labels = c('Missense', 'Truncating', 'In_Frame_InDel')),
                show_column_names = TRUE,
                show_row_names = T,
                pct_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12)) 

png('./Results/tmp1.png',width = 5.4, height = 2.5, res = 300, units = 'in')
p1
dev.off()
png('./Results/tmp2.png',width = 8, height = 2.5, res = 300, units = 'in')
p2
dev.off()
png('./Results/tmp3.png',width = 4.5, height = 2.5, res = 300, units = 'in')
p3
dev.off()

img1 <- image_read('./Results/tmp1.png')
img2 <- image_read('./Results/tmp2.png')
img3 <- image_read('./Results/tmp3.png')
file.remove('./Results/tmp1.png','./Results/tmp2.png','./Results/tmp3.png')

combined <- image_append(c(img1, img2, img3))
image_write(combined, './Results/Figure 2D.png')



