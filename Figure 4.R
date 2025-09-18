rm(list = ls())
setwd(getwd())
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(Seurat)
library(cytosignal)
library(Matrix)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 4A
plot_mat <- readRDS('./Data/Figure 4A.rds')

p <- ggplot(plot_mat, aes(x = log10p.adjust, y = Description, fill = log10p.adjust)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_fill_gradientn(colours = c(brewer.pal(7, 'Set1')[6], brewer.pal(8, 'Set1')[1])) +
  theme_bw()+
  # theme_void()+
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = 'black'), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6))

ggsave('./Results/Figure 4A.png', p, width = 4, height = 3)


# Figure 4B
load('./Data/Figure 4B.rda')

p <- ggplot(MPs_corr, aes(x = MP_type1, y = MP_type2, size = spearman_p.adjust_inv, color =  spearman_coef)) + 
  scale_color_gradientn(colours = c(colorRampPalette(c('#2171b5', '#ffffbf'))(709), 
                                    colorRampPalette(c('#ffffbf', '#d73027'))(300), 
                                    colorRampPalette(c('#d73027', '#bd0026'))(221))) + 
  geom_point(stat = 'identity') + 
  xlab('') + ylab('') + 
  ggtitle('Spearman correlation of MPs') + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./Results/Figure 4B.MPs_correlation.png', p, width = 5.65, height = 4)

png('./Results/Figure 4B.MPs_pathology_Roe.png', width = 7, height = 3.5, res = 300, units = 'in')
pheatmap(MPs_pathol_Roe, color = c(colorRampPalette(c('#4575b4', 'white'))(100), colorRampPalette(c('white', '#d73027'))(250)), 
         cluster_rows = F, cluster_cols = F)
dev.off()

png('./Results/Figure 4B.MPs_lesion_type_Roe.png', width = 5.5, height = 1.8, res = 300, units = 'in')
pheatmap(MPs_lesion_Roe, color = c(colorRampPalette(c('#313695', 'white'))(760), colorRampPalette(c('white', 'orange'))(481)), 
             cluster_rows = F, cluster_cols = F)
dev.off()


# Figure 4C
load('./Data/Figure 4C.rda')

p1 <- ggplot(barplot_mat, aes(x = Pairs, y = Freq, fill = Group) ) +
  geom_bar(width = 0.7, stat = 'identity') +
  scale_fill_manual(values = c('#4DBBD5FF', '#11A579', '#E64B35FF')) +
  labs(x = '', y = 'Sample Number') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'right') 

p2 <- ggplot(piechart_mat, aes(x = '', y = percentage, fill = Group) ) + 
  scale_fill_manual(values = c('#11A579', '#66C5CC', '#F2B701', '#E73F74')) + 
  geom_bar(width = 1, stat = 'identity') +
  coord_polar('y', start = 0) +
  labs(x = '', y = '') +
  theme_void() +
  theme(axis.text.x = element_text(angle = 0, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'right')

png('./Results/Figure 4C.png', width = 5, height = 1.7, res = 300, units = 'in')
plot_grid(p1, p2, align = 'ht', ncol = 2)
dev.off()


# Figure 4D
seurat_obj <- readRDS('./Data/Figure 4D.rds')

png('./Results/Figure 4D.png', width = 3.5, height = 3, res = 300, units = 'in')
VlnPlot(object = seurat_obj, group.by = 'celltype', features = 'IL1R1', pt.size = F, col = c('AT1' = '#1f77b4', 'AT2' = '#aec7e8', 'AIC' = '#e377c2', 'KAC' = '#d62728', 'Tumor' = '#f7b6d2'))
dev.off()


# Figure 4E
load('./Data/Figure 4E.rda')

intr.use <- names(pair)

p1 <- plotEdge(P24_AAH_cs, intr.use, slot.use = 'GauEps-Raw', pt.size = 0.25, colors.list = P24_AAH_Lable_colors) 
p2 <- plotEdge(P24_LUAD_cs, intr.use, slot.use = 'GauEps-Raw', pt.size = 0.25, colors.list = P24_LUAD_Lable_colors) 

png('./Results/Figure 4E.AAH.png', width = 8, height = 7, res = 300, units = 'in')
print(p1) 
dev.off()

png('./Results/Figure 4E.LUAD.png', width = 8, height = 7, res = 300, units = 'in')
print(p2)
dev.off()


# Figure 4G was created by using 'Figure 4G.py'
