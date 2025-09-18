rm(list = ls())
setwd(getwd())
library(Seurat)
library(ggplot2)
library(monocle3)
library(ggsci)
library(viridis)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggalluvial)

if(!dir.exists('./Results')){
    dir.create('./Results')
}

# Figure 3A
Epi <- readRDS('./Data/Figure 3A.rds')

png('./Results/Figure 3A.png', width = 5, height = 5, res = 300, units = 'in')
DimPlot(Epi, reduction = 'umap', label = TRUE, group.by = 'celltype', pt.size = 1.5) + scale_color_manual(values = c(Basal = '#1f77b4', AT1 = '#ff7f0e', AT2 = '#e377c2', KAC = '#2ca02c', Ciliated = '#ffbb78', Club = '#ff9896', Club.secretory = '#9467bd', AIC = '#d62728', Tumor = '#f7b6d2'))
dev.off()


# Figure 3B-D
monocle3_obj<-readRDS('./Data/Figure 3B-D.downsample.rds')

p1 <- plot_cells(monocle3_obj, color_cells_by = 'celltype', cell_size = 0.6, show_trajectory_graph = FALSE, group_label_size = 0.5, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_groups_by_cluster = F) + ggtitle(NULL) + scale_color_manual(values = c(AT1 = '#aec7e8', AT2 = '#e377c2', AIC = '#2ca02c', KAC = '#d62728', 'Precursor' = '#bcbd22', 'Invasive' = '#17becf'))
p1$layers[[1]]$aes_params$colour <- 'transparent'
p2 <- plot_cells(monocle3_obj, color_cells_by = 'Pseudotime', cell_size = 0.6, show_trajectory_graph = TRUE, trajectory_graph_segment_size = 0.25, group_label_size = 0.5, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p2$layers[[1]]$aes_params$colour <- 'transparent'
p3 <- plot_cells(monocle3_obj, color_cells_by  = 'CytoTRACE', cell_size = 0.6, show_trajectory_graph = FALSE, group_label_size = 0.5, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + ggtitle(NULL) + scale_color_gradientn(colors = turbo(100))
p3$layers[[1]]$aes_params$colour <- 'transparent'

png(paste0('./Results/Figure 3B-D.png'), width = 18, height = 6, res = 300, units = 'in')
plot_grid(p1, p2, p3, align = 'hv', ncol = 3)
dev.off()


# Figure 3E
load('./Data/Figure 3E.rda')
custom_magma <- c(colorRampPalette(c('white', rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

png('./Results/Figure 3E.png', width = 10, height = 9, units = 'in', res = 300)
pheatmap(MP_mat, color = custom_magma, breaks = seq(4, 30, by = 0.1), 
         scale = 'none', cluster_rows = F, cluster_cols = F, 
         annotation_colors = column_annotation_colors, 
         annotation_col = column_annotation, border_color = NA, 
         show_rownames = F, show_colnames = F, 
         fontsize_row = 3)
dev.off() 


# Figure 3F
plot_mat <- readRDS('./Data/Figure 3F.rds')

p <- ggplot(plot_mat, aes(x = '', y = percentage, fill = MP)) + 
  geom_bar(width = 1, stat = 'identity') +
  coord_polar('y', start = 0) +
  scale_fill_manual(values = c(MP1 = '#F6CF71', MP2 = '#3969AC', MP3 = '#80BA5A', MP4 = '#F2B701', MP5 = '#11A579', MP6 = '#CF1C90', MP7 = '#66C5CC', MP8 = '#f97b72', MP9 = '#ed5887')) + 
  facet_wrap(vars(CellType), ncol = 4) +
  labs(x = '', y = '') +
  theme_void()  +
  theme(axis.text.x = element_text(angle = 0, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'right') 

ggsave('./Results/Figure 3F.png', p, width = 10, height = 5)


# Figure 3G
MP_mat<-readRDS('./Data/Figure 3G.rds')
custom_magma <- c(colorRampPalette(c('white', rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

png('./Results/Figure 3G.png', width = 3.3, height = 2, res = 300, units = 'in')
pheatmap(MP_mat, color = custom_magma, breaks = seq(4, 25, by = 0.1), 
         scale = 'none', 
         cluster_cols = F, cluster_rows = F, 
         show_rownames = T, show_colnames = F, 
         fontsize_row = 5, fontsize_col = 5, 
         border_color = NA)
dev.off()


# Figure 3H
plot_mat <- readRDS('./Data/Figure 3H.rds')

p <- ggplot(plot_mat, aes(x = Group_ID, y = percentage, fill = Group_ID)) + 
  scale_fill_manual(values = c(Normal = '#11A579', AAH = '#66C5CC', AIS = '#F6CF71', MIA = '#f97b72', LUAD = '#ed5887')) + 
  geom_boxplot(lwd = 0.1, width = 0.8, show.legend = F, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.1) +
  facet_wrap(vars(MP), ncol = 5, scales = 'free') +
  stat_compare_means(label = 'p.format', vjust = 0, size = 2)+  
  labs(x = '', y = 'MP% in each sample of snRNA-Seq data') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none') 

ggsave('./Results/Figure 3H.png', p, width = 5, height = 3)


# Figure 3I
plot_mat <- readRDS('./Data/Figure 3I.rds')

p <- ggplot(plot_mat, aes(x = Patient_ID, y = ifelse(Malignant_State == 'Malignant', -Freq, Freq), fill = MP)) +
  geom_bar(stat = 'identity') +                            
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = c(MP1 = '#F6CF71', MP2 = '#3969AC', MP3 = '#80BA5A', MP4 = '#F2B701', MP5 = '#11A579', MP6 = '#CF1C90', MP7 = '#66C5CC', MP8 = '#f97b72', MP9 = '#ed5887')) + 
  labs(x = '', y = 'Invasive (bottom)  ---  Non-invasive (top)') +
  theme_classic()+
  geom_hline(yintercept = 0, color = 'white')

ggsave('./Results/Figure 3I.png', p, width = 8, height = 5)


# Figure 3J
plot_mat <- readRDS('./Data/Figure 3J.rds')

p <- ggplot(plot_mat, aes(x = MP_label, y = percentage, fill = MP_label)) + 
  scale_fill_manual(values = c(MP1 = '#F6CF71', MP2 = '#3969AC', MP3 = '#80BA5A', MP4 = '#F2B701', MP5 = '#11A579', MP6 = '#CF1C90', MP7 = '#66C5CC', MP8 = '#f97b72', MP9 = '#ed5887')) + 
  geom_boxplot(lwd = 0.2, width = 0.8, show.legend = F, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.1) +
  labs(x = '', y = 'PRII MPs (%) by sample') +
  stat_compare_means(label = 'p.format', vjust = 0, size = 2)+  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black', ), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none')

ggsave('./Results/Figure 3J.png', p, width = 2.1, height = 2)


# Figure 3K were created by using 'Figure 3K 3N.py'

# Figure 3L
plot_mat <- readRDS('./Data/Figure 3L.rds')

p <- ggplot(plot_mat, aes(x = Pathology, stratum = MP, alluvium = MP, y = percentage, fill = MP)) +
  geom_stratum(width = 0.5, col=NA) +
  geom_alluvium(width = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c(MP1 = '#F6CF71', MP2 = '#3969AC', MP3 = '#80BA5A', MP4 = '#F2B701', MP5 = '#11A579', MP6 = '#CF1C90', MP7 = '#66C5CC', MP8 = '#f97b72', MP9 = '#ed5887')) + 
  labs(x = '', y = 'Fraction of each MP (%)')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, size = 12,  color = 'black', hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 12,  color = 'black'))

ggsave('./Results/Figure 3L.png', p, width = 3.33, height = 3)


# Figure 3M
load('./Data/Figure 3M.rda')

p1 <- ggplot(RP_clone_mat, aes(x = clone_new, y = percentage, fill = lesions)) + 
  scale_fill_manual(values = c('#11A579', '#66C5CC', '#ed5887')) + 
  geom_boxplot(lwd = 0.1, width = 0.7, show.legend = T, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(0.8), dotsize = 0.8)+
  ggtitle('Evolution Pattern 1a') +
  labs(x = '', y = 'KAC score') +
  stat_compare_means(comparisons = list(c('Shared', 'Invasive-specific')), method = 't.test', label = 'p.format', vjust = 0, size = 4)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black', ), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'right') 

p2 <- ggplot(KAC_mat, aes(x = clone_new, y = KAC_Score, fill = clone_new)) +
  geom_violin(lwd = 0.2) +
  geom_boxplot(lwd = 0.2, width = 0.3, fill = 'white', outlier.shape = NA) +
  scale_fill_manual(values = c('#1f77b4', '#e377c2', '#aec7e8', '#d62728')) +
  stat_compare_means(comparisons = list(c('Ref', 'Shared'), c('Ref', 'Invasive-specific'), c('Shared', 'Invasive-specific')), method = 't.test', label = 'p.format', vjust = 0, size = 1.5) +
  ggtitle('Evolution Pattern 1a') +
  labs(x = '', y = 'RPII clones (%) by sample') +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 20, size = 6, color = 'black', vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        strip.text.x = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.text.y = element_text(size = 6, color = 'black', face = 'bold'), 
        strip.background = element_rect(fill = NA, colour = NA), 
        text = element_text(size = 6), legend.position = 'none')

png('./Results/Figure 3M.png', width = 7, height = 3.5, res = 300, units = 'in')
plot_grid(p1, p2, align = 'hv', ncol = 2)
dev.off()


# Figure 3N were created by using 'Figure 3K 3N.py'