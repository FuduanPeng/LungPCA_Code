rm(list = ls())
setwd(getwd())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(cowplot)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 6B, 6D
seurat_obj <- readRDS('./Data/Figure 6B 6D.rds')

# Figure 6B
png('./Results/Figure 6B.png', width = 5.5, height = 4, res = 300, units = 'in')
DimPlot(seurat_obj, reduction = 'umap', label = F, group.by = 'celltype', pt.size = 0.1) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club and Secretory' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'Proliferating' = '#ff7f0e', 'Tuft' = '#1f77b4', 'Tumor' = '#f7b6d2'))
dev.off()

# Figure 6D
Epi_alveolar <- seurat_obj[, seurat_obj$celltype %in% c('KAC', 'AT1', 'AT2', 'Tumor')]
png('./Results/Figure 6D.png', width = 5, height = 6, res = 300, units = 'in')
VlnPlot(object = Epi_alveolar, features = c('Il1r1', 'Rela', 'Relb', 'Nfkb1'), group.by = 'celltype', pt.size = F, ncol = 2, col = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Tumor' = '#f7b6d2'))
dev.off()


# Figure 6C
plot_mat <- readRDS('./Data/Figure 6C.rds')

p <- ggplot(plot_mat, aes(x = log10p.adjust, y = Description, fill = log10p.adjust)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_fill_gradientn(colours = c(brewer.pal(7, 'Set1')[6], brewer.pal(8, 'Set1')[1])) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = 'black'), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6))

ggsave('./Results/Figure 6C.png', p, width = 4.2, height = 3)


# Figure 6E and 6F were created by using 'Figure 6E 6F 7D.py'

# Figure 6G were created by using 'Figure 6G.py'


# Figure 6J, 6K
load('./Data/Figure 6J 6K.rda')

p1 <- ggplot(number_mat, aes(x = Group_ID, y = Node_Number, fill = Group_ID)) + 
  geom_boxplot(lwd = 0.2, width = 0.5, show.legend = F, outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.5), size = 0.5) +
  scale_fill_manual(values = c('#17becf', '#ff7f0e', '#e377c2')) +
  labs(x = '', y = 'Organoid number') +
  stat_compare_means(comparisons = list(c('Control', 'IL-1β' ), c('Control', 'IM' ), c('IL-1β', 'IM' )), method = 't.test', label = 'p.format', vjust = 0, size = 4) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = 'black', vjust = 0.95), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        strip.text.x = element_text(size = 12, color = 'black', face = 'bold'), 
        strip.text.y = element_text(size = 12, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        legend.position = 'none') 

p2 <- ggplot(size_mat, aes(x = Group_ID, y = Node_Size, fill = Group_ID)) + 
  geom_boxplot(lwd = 0.2, width = 0.5, show.legend = F, outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.3), size = 0.5) +
  scale_fill_manual(values = c('#17becf', '#ff7f0e', '#e377c2')) +
  labs(x = '', y = 'Organoid size (um)') +
  stat_compare_means(comparisons = list(c('Control', 'IL-1β' ), c('Control', 'IM' ), c('IL-1β', 'IM' )), method = 't.test', label = 'p.format', vjust = 0, size = 4) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = 'black', vjust = 0.95), 
        axis.text.y = element_text(size = 12, color = 'black'), 
        strip.text.x = element_text(size = 12, color = 'black', face = 'bold'), 
        strip.text.y = element_text(size = 12, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        legend.position = 'none') 

png('./Results/Figure 6J 6K.png', width = 3*2, height = 4, res = 300, units = 'in')
plot_grid(p1, p2, align = 'hv', nrow = 1)
dev.off()


