rm(list = ls())
setwd(getwd())
library(ggplot2)
library(cowplot)
library(ggsankey)
library(dplyr)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 1B
load('./Data/Figure 1B.rda')

p1<-ggplot(sankeypt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) + 
  geom_sankey(flow.alpha = 0.5, node.color = 1) +  
  geom_sankey_label(size = 3, color = 1.2) + 
  scale_fill_manual(values = c(Epi_spots = '#1f78b4', nonEpi_spots = '#a6cee3', Normal = '#33a02c', AAH = '#b2df8a', AIS = '#fdbf6f', MIA = '#fb9a99', LUAD = '#ff7f00', 'Tumor necrosis' = '#e31a1c', Lymphoid = '#ffff99', Myeloid = '#6a3d9a', Fibroblast = '#cab2d6', Vessel = '#b15928', 'Smooth muscle cells' = '#a6761d', Other = '#d9d9d9')) +
  theme_sankey(base_size = 6) + 
  theme(legend.position = 'none') + 
  xlab('') 

p2<-ggplot(data = barpt, aes(x = Spots_Number, y = Spots_Type, fill = Spots_Type)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = c(Normal = '#33a02c', AAH = '#b2df8a', AIS = '#fdbf6f', MIA = '#fb9a99', LUAD = '#ff7f00', 'Tumor necrosis' = '#e31a1c', Lymphoid = '#ffff99', Myeloid = '#6a3d9a', Fibroblast = '#cab2d6', Vessel = '#b15928', 'Smooth muscle cells' = '#a6761d', Other = '#d9d9d9')) +
  ylab('') + 
  xlab('log10(Spots Number)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        strip.text.x = element_text(size = 6), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), 
        legend.position = 'none') 

p3<-ggplot(data = boxpt, aes(x = Spots_Number, y = Spots_Type, col = Spots_Type, fill = Spots_Type)) + 
  geom_boxplot(lwd = 0.2, alpha = 0.8, width = 0.5, show.legend = F, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1)+
  scale_color_manual(values = c(Normal = '#33a02c', AAH = '#b2df8a', AIS = '#fdbf6f', MIA = '#fb9a99', LUAD = '#ff7f00', 'Tumor necrosis' = '#e31a1c', Lymphoid = '#ffff99', Myeloid = '#6a3d9a', Fibroblast = '#cab2d6', Vessel = '#b15928', 'Smooth muscle cells' = '#a6761d', Other = '#d9d9d9')) +
  scale_fill_manual(values = c(Normal = '#33a02c', AAH = '#b2df8a', AIS = '#fdbf6f', MIA = '#fb9a99', LUAD = '#ff7f00', 'Tumor necrosis' = '#e31a1c', Lymphoid = '#ffff99', Myeloid = '#6a3d9a', Fibroblast = '#cab2d6', Vessel = '#b15928', 'Smooth muscle cells' = '#a6761d', Other = '#d9d9d9')) +
  ylab('') + 
  xlab('log10(Spots Number)\nper sample') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_blank(), 
        strip.text.x = element_text(size = 6), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), 
        legend.position = 'none') 

png('./Results/Figure 1B.png', width = 4*3, height = 5, res = 300, units = 'in')
plot_grid(p1, p2, p3, ncol = 3, 
          rel_widths = c(2, 0.8, 0.8), 
          align = 'hv')
dev.off()


# Figure 1C
seurat_obj<-readRDS('./Data/Figure 1C.downsample_25_percent.rds')

p1<-DimPlot(dseurat_obj, reduction = 'umap', label = TRUE, group.by = 'Pathology', pt.size = 0.8) + scale_color_manual(values = c('#1F78B4', '#FF7F00', '#FB9A99', '#B2DF8A'))
p2<-FeaturePlot(seurat_obj, features = c('KRT8', 'SFTPC', 'MUC5B', 'EGFR', 'CEACAM5', 'TFF3', 'CLDN4', 'KRT7', 'CEACAM6'), ncol = 3) & scale_colour_gradientn(colours = c('white', '#ffffe5', '#ffffcc', '#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026'))

png('./Results/Figure 1C.png', width = 7.5*2, height = 6, res = 300, units = 'in')
plot_grid(p1, p2, align = 'hv', ncol = 2)
dev.off()


# Figure 1D
scaled_data<-readRDS('./Data/Figure 1D.rds')

png('./Results/Figure 1D.png', width = 11, height = 2, res = 300, units = 'in')
pheatmap(scaled_data, cluster_rows = F, cluster_cols = F, fontsize = 6, border_color = 'white', angle_col = 45, color = colorRampPalette(rev(brewer.pal(n = 10, name  = 'RdBu')))(256))
dev.off()


# Figure 1E was created by using 'Figure 1E 4F.py'

