rm(list = ls())
setwd(getwd())
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 5A 
seurat_obj_downsample <- readRDS('./Data/Figure 5A.rds')

png('./Results/Figure 5A.downsampled.png', width = 10, height = 8, res = 300, units = 'in')
DimPlot(seurat_obj_downsample, reduction = 'umap_harmony', group.by = 'Major_cell_type', label = F, raster = T, 
        cols = c('#FF0000', '#FF00FF', '#808000',  '#800000', '#54278f', '#FF9900', '#33FFFF', '#E6BEFF', '#66FF33', 
               '#A9A9A9', '#FABEBE',  '#008080', '#b3de69', '#FFFF00', '#8dd3c7', '#D2691E', '#4363D8', '#deebf7'))
dev.off()


# Figure 5B-C
seurat_obj <- readRDS('./Data/Figure 5B 5C.rds')

# Figure 5B
png('./Results/Figure 5B.png', width = 10, height = 8, res = 300, units = 'in')
DimPlot(seurat_obj, reduction = 'umap_harmony', group.by = 'lineage', label = F, raster = T, 
        cols = c('#0570b0', '#8dd3c7', '#006837', '#feb24c', '#bd0026', '#4eb3d3', 
               '#b3de69', '#fcc5c0', '#54278f', '#fddbc7', '#fccde5'))
dev.off()

# Figure 5C
p1<-VlnPlot(seurat_obj, features = 'KAC_signature', group.by = 'lineage', slot = 'data', pt.size = 0, 
        cols = c('#0570b0', '#8dd3c7', '#006837', '#feb24c', '#bd0026')) + 
  stat_summary(fun = median, geom = 'point', shape = 21, color = 'black', size = 3, fill = 'black')

p2<-VlnPlot(seurat_obj, features = 'Inflammatory_pathway', group.by = 'lineage', slot = 'data', pt.size = 0, 
        cols = c('#0570b0', '#8dd3c7', '#006837', '#feb24c', '#bd0026'))+
    stat_summary(fun = median, geom = 'point', shape = 21, color = 'black', size = 3, fill = 'black')

png('./Results/Figure 5C.png', width = 4, height = 6, res = 300, units = 'in')
p1+p2
dev.off()


# Figure 5D
plot_mat <- readRDS('./Data/Figure 5D.rds')
p <- ggplot(plot_mat, aes(center_cell_type, neighborhood_cell_prop, fill = neighborhood_cell_type)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = c('#FF0000', '#FF00FF', '#808000', '#800000', '#54278f', '#FF9900', '#33FFFF', '#E6BEFF', '#66FF33', 
                             '#A9A9A9', '#FABEBE', '#008080', '#b3de69', '#FFFF00', '#8dd3c7', '#D2691E', '#4363D8', '#deebf7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, color = 'black', hjust = 1), 
        axis.text.y = element_text(color = 'black'), 
        axis.ticks.length = unit(0.1, 'cm'), 
        legend.position = 'right') +
  guides(fill = guide_legend(title = NULL))
ggsave('./Results/Figure 5D.png', p, width = 5, height = 5)


# Figure 5E
plot_mat <- readRDS('./Data/Figure 5E.rds')

p <- ggplot(plot_mat, aes(center_cell_type, Myeloid_C15, fill = center_cell_type)) +
  geom_violin(scale = 'width', adjust = 15) +
  geom_boxplot(width = 0.3, fill = 'white', outlier.alpha = 0) +
  scale_fill_manual(values = c('#0570b0', '#fb8072', '#006837', '#feb24c', '#4eb3d3', '#bd0026', '#b3de69', '#8dd3c7', '#fcc5c0', 
                             '#fddbc7', '#9A6324', '#E6194B', '#3CB44B', '#4363D8', '#e21c4d', '#3eb34b', '#4f66ae')) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, color = 'black', hjust = 1), 
        axis.text.y = element_text(color = 'black'), 
        axis.ticks.length = unit(0.1, 'cm'), 
        legend.position = 'right') +
  guides(fill = guide_legend(title = NULL))
ggsave('./Results/Figure 5E.png', width = 5, height = 5)


### Figure 5F was created using Xenium Explorer


# Figure 5G
block_info <- readRDS('./Data/Figure 5G.rds')

plot_mat <- as.data.frame(table(block_info$Pathology))
colnames(plot_mat) <- c('Pathology', 'Num')
plot_mat$Pathology <- factor(plot_mat$Pathology, levels = c('No_lesion', 'AAH', 'AIS', 'MIA', 'LUAD_lepidic', 'LUAD_acinar', 'LUAD_papillary', 'LUAD_solid', 'Unassigned'))

p <- ggplot(plot_mat, aes(x = '', y = Num, fill = Pathology)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar('y', start = 0)+
  scale_fill_manual(values = c('#B0C4DE', '#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#46F0F0', '#d9d9d9')) +
  theme_classic()
ggsave('./Results/Figure 5G.png', p, width = 5, height = 5)


# Figure 5H
seurat_obj <- readRDS('./Data/Figure 5H.rds')

png('./Results/Figure 5H.png', width = 10, height = 8, res = 300, units = 'in')
DimPlot(seurat_obj, reduction = 'umap_harmony', label = F, raster = T, group.by = 'Major_cell_type', 
        cols = c('#FF0000', '#FF00FF', '#808000',  '#800000', '#54278f', '#FF9900', '#33FFFF', '#66FF33',  
               '#FABEBE',  '#008080', '#b3de69', '#FFFF00', '#8dd3c7', '#4363D8', '#deebf7'))
dev.off()


# Figure 5I was created using Xenium Explorer


# Figure 5J
plot_mat <- readRDS('./Data/Figure 5J.rds')

p <- ggplot(plot_mat, aes(center_cell_type, neighborhood_cell_prop, fill = neighborhood_cell_type)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = c('#FF0000', '#FF00FF', '#808000',  '#800000', '#54278f', '#FF9900', '#33FFFF', '#66FF33', 
                             '#FABEBE',  '#008080', '#b3de69', '#FFFF00', '#8dd3c7', '#4363D8', '#deebf7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, color = 'black', hjust = 1), 
        axis.text.y = element_text(color = 'black'), 
        axis.ticks.length = unit(0.1, 'cm'), 
        legend.position = 'right') +
  guides(fill = guide_legend(title = NULL))
ggsave('./Results/Figure 5J.png', p, width = 5, height = 5)


### Figure 5K
plot_mat <- readRDS('./Data/Figure 5K.rds')

p <- ggplot(plot_mat, aes(Pathology2, Mac_C4_IL1B, fill = Pathology2)) +
  geom_violin(scale = 'width', adjust = 15) +
  geom_boxplot(width = 0.3, fill = 'white', outlier.alpha = 0) +
  scale_fill_manual(values = c('#E6194B', '#3CB44B', '#FFE119', '#4363D8')) +
  theme_classic() +                            
  theme(axis.text.x = element_text(angle = 45, color = 'black', hjust = 1), 
        axis.text.y = element_text(color = 'black'), 
        axis.ticks.length = unit(0.1, 'cm'), 
        legend.position = 'right') +
  guides(fill = guide_legend(title = NULL))
ggsave('./Results/Figure 5K.png', p, width = 4, height = 4)


# Figure 5L was created using Xenium Explorer

