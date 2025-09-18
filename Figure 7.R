rm(list = ls())
setwd(getwd())
library(Seurat)
library(cowplot)
library(ggpubr)
library(dplyr)
library(ggplot2)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 7B and 7C were created using Prism

# Figure 7D were created using 'Figure 6E 6F 7D.py'

# Figure 7E, 7F
load('./Data/Figure 7E 7F 7G.rda')

# Figure 7E
png('./Results/Figure 7E.png', width = 5.5, height = 4.2, res = 300, units = 'in')
DimPlot(seurat_obj, reduction = 'umap', label = T, group.by = 'celltype', pt.size = 0.05) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'AIC' = '#1f77b4', 'Tumor cells' = '#f7b6d2')) + ggtitle('')
dev.off()

# Figure 7F
p1 <- ggplot(mat_3mo, aes(x = Group_ID, y = percentage, fill = Group_ID)) + 
  geom_boxplot(lwd = 0.2, alpha = 1, width = 0.7,  show.legend = F, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5)+
  scale_fill_manual(values = c('#66C5CC', '#F2B701', '#7F3C8D', '#CF1C90')) + 
  ggtitle('3 months') +
  labs(x = '', y = 'Fraction in Epithelial cells') +
  stat_compare_means(comparisons = list(c('anti-IL-1β', 'IgG'), c('anti-PD-1', 'IgG'), c('Combo', 'IgG'), c('anti-IL-1β', 'Combo'), c('Combo', 'anti-PD-1')), method = 't.test', label = 'p.format', vjust = 0, size = 1.5)+  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none') 

p2 <- ggplot(mat_7mo, aes(x = Group_ID, y = percentage, fill = Group_ID)) + 
  geom_boxplot(lwd = 0.2, alpha = 1, width = 0.7,  show.legend = F, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5)+
  scale_fill_manual(values = c('#66C5CC', '#F2B701', '#7F3C8D', '#CF1C90')) + 
  ggtitle('3 months') +
  labs(x = '', y = 'Fraction in Epithelial cells') +
  stat_compare_means(comparisons = list(c('anti-IL-1β', 'IgG'), c('anti-PD-1', 'IgG'), c('Combo', 'IgG'), c('anti-IL-1β', 'Combo'), c('Combo', 'anti-PD-1')), method = 't.test', label = 'p.format', vjust = 0, size = 1.5)+  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none') 

png('./Results/Figure 7F.png', width = 4, height = 3.35, res = 300, units = 'in')
plot_grid(p1, p2, align = 'hv', nrow = 1)
dev.off()

# Figure 7G
seurat_obj <- seurat_obj[, seurat_obj$celltype %in% c('AT1', 'AT2', 'AIC', 'KAC', 'Tumor cells')]
p1 <- DimPlot(seurat_obj[, seurat_obj$treatment == 'IgG'], reduction = 'umap', label = F, group.by = 'celltype', pt.size = 0.01) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'AIC' = '#1f77b4', 'Tumor cells' = '#f7b6d2')) + ggtitle('IgG')
p2 <- DimPlot(seurat_obj[, seurat_obj$treatment == 'anti-IL-1β'], reduction = 'umap', label = F, group.by = 'celltype', pt.size = 0.01) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'AIC' = '#1f77b4', 'Tumor cells' = '#f7b6d2')) + ggtitle('anti-IL-1β')
p3 <- DimPlot(seurat_obj[, seurat_obj$treatment == 'anti-PD-1'], reduction = 'umap', label = F, group.by = 'celltype', pt.size = 0.01) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'AIC' = '#1f77b4', 'Tumor cells' = '#f7b6d2')) + ggtitle('anti-PD-1')
p4 <- DimPlot(seurat_obj[, seurat_obj$treatment == 'Combo'], reduction = 'umap', label = F, group.by = 'celltype', pt.size = 0.01) + scale_color_manual(values = c('KAC' = '#17becf', 'AT1' = '#aec7e8', 'AT2' = '#e377c2', 'Basal' = '#98df8a', 'Ciliated' = '#ffbb78', 'Club' = '#ff9896', 'Neuroendocrine' = '#9467bd', 'AIC' = '#1f77b4', 'Tumor cells' = '#f7b6d2')) + ggtitle('Combo')

png('./Results/Figure 7G.png', width = 5.5*2, height = 4*2, res = 300, units = 'in')
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()


# Figure 7H was created using ImageScope


# Figure 7I
plot_mat <- readRDS('./Data/Figure 7I.rds')

p <- ggplot(plot_mat, aes(x = Group, y = macro.Frac, fill = Group)) +
  geom_boxplot(lwd = 0.1, width = 0.7, show.legend = T, outlier.shape = NA, position = position_dodge(0.8)) +
  scale_fill_manual(values = c('#66C5CC', '#F2B701', '#7F3C8D', '#CF1C90')) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5) +
  labs(x = '', y = 'Macrophage Fraction') +
  stat_compare_means(comparisons = list(c('IgG', 'anti-IL-1β'), c('IgG', 'anti-PD-1'), c('IgG', 'Combo'), c('Combo', 'anti-PD-1')), method = 't.test', label = 'p.format', vjust = 0, size = 1.5)+  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none') 

ggsave('./Results/Figure 7I.png', p, width = 1.3, height = 1.75)






