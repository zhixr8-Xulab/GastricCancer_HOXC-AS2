#==========Load required packages==============
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggprism)
library(scales)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggridges)

meta <- readRDS(file = 'meta.rds')

meta$celltype <- factor(meta$celltype, levels = c('Epithelials', 'Endocrines', 'Ehdothelials', 'Fibroblasts', 'Pericytes', 'Mesothelials',
                                                  'Mast cells', 'Macro/Mono', 'T cells', 'B cells', 'Plasmas'))

meta$epitype <- factor(meta$epitype, levels = c('Epi_01', 'Epi_02', 'Epi_03', 'Epi_04', 'Epi_05',
                                                'Epi_06', 'Epi_07', 'Epi_08', 'Epi_09', 'Epi_10'))

meta$fibrotype <- factor(meta$fibrotype, levels = c('Fibro_C1', 'Fibro_C2', 'Fibro_C3', 'Fibro_C4', 'Fibro_C5', 'Fibro_C6', 'Fibro_C7',
                                                    'Peri_C1', 'Peri_C2', 'Peri_C3', 'Mesothelial'))

celltype_col <- c("#00A087B2", "#91D1C2B2", "#8491B4B2", '#e3843b', "#3C5488B2", "#4DAF4A",
                  "#4DBBD5B2", "#F39B7FB2", "#E64B35B2", "#CC79A7", "#A14462")

tissue_col <- c('#1b5dbd', '#789dcb', '#f46891', '#c31e3a')
  
fibrotype_col <- c('#f29145', '#fec888', '#e2a09c', '#c2b2d0', '#c54a85', '#be3f33',
                   '#cf91bd','#2171b5', '#b3cbde', '#6baed6', '#c0db9a')

epitype_col <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                 "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3")


#==========Fig 3f============
pdf('fig_3f.pdf', height = 5, width = 4.5)
meta %>% filter(celltype %in% c('Fibroblasts', 'Pericytes', 'Mesothelials')) %>%
  ggplot(., aes(x = tissue, y = AGCRS, fill = tissue)) + 
  geom_boxplot(outlier.alpha = 0) +
  facet_wrap(~celltype, nrow = 1) +
  theme_prism(border = T, base_line_size = 0.4, base_size = 12, axis_text_angle = 45) +
  guides(fill = 'none') + labs(x = '', y = 'AGCRS Score') +
  scale_fill_manual(values = tissue_col)
dev.off()

#==========Fig 3g==========
dat <- meta %>% filter(tissue == 'Primary') %>% filter(celltype %in% c('Fibroblasts', 'Pericytes', 'Mesothelials')) %>%
  filter(stage %in% c('I', 'II', 'III', 'IV'))
pdf('fig_3g.pdf', height = 5.5, width = 5.5)
ggplot(dat, aes(x = stage, y = AGCRS)) +
  geom_violin(aes(fill = stage), position = position_dodge(width = 1)) +
  stat_summary(aes(fill = stage), fun = median, position = position_dodge(width = 1), size = 0.3) +
  labs(x = 'Stage', y = 'Primary tissue AGCRS score') +
  guides(fill = 'none') +
  scale_y_continuous(breaks = seq(0,1,.5)) +
  facet_grid(~ celltype, scales = 'free_x', space = 'free_x') +
  # scale_fill_manual(values = c('#e3843b', '#96a9d3', '#519d52')) +
  scale_fill_manual(values = rev(pal_npg()(4))) +
  theme_prism(base_line_size = 0.6) +
  theme(legend.position = 'top', strip.background = element_rect(fill = NA))
dev.off()

#==========fig 3h==========
fibro_meta <- meta %>% filter(celltype %in% c('Fibroblasts', 'Pericytes', 'Mesothelials'))
fibro_prop <- table(fibro_meta$fibrotype, fibro_meta$tissue) %>% prop.table(., 2) %>% as.data.frame()

pdf('fig_3h.pdf', height = 8, width = 6.5)
fibro_prop %>% filter(Var1 %in% c('Fibro_C1', 'Fibro_C2', 'Fibro_C3', 'Fibro_C4', 'Fibro_C5', 'Fibro_C6', 'Fibro_C7')) %>%
  ggplot(., aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = 'identity', position = 'fill', width = 0.7, color = 'black') +
  labs(x = '', y = 'Percentage(%)') +
  scale_fill_manual(values = tissue_col) +
  theme_prism(border = T, base_line_size = 0.4, base_fontface = 'bold') +
  scale_y_continuous(expand = c(0.01,0.01), labels = seq(0,100,25)) +
  theme(legend.position = 'right', legend.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

#==========Fig 3i==========
dat <- meta %>% filter(celltype == 'Fibroblasts')
dat$AGCRS <- rescale(dat$AGCRS, to = c(0, 1))
pdf('fig_3i.pdf', height = 5, width = 9)
ggplot(dat, aes(x = tissue, y = AGCRS, fill = tissue)) + 
  geom_boxplot(outlier.alpha = 0, position = position_dodge(0.8), width = 0.5) +
  facet_grid(~fibrotype, scales = 'free_x', space = 'free_x') +
  guides(fill = 'none') + labs(x = '', y = 'AGCRS Score') +
  scale_fill_manual(values = tissue_col) +
  theme_prism(border = T, base_size = 12, axis_text_angle = 45)
dev.off()

#===========Fig 4h=========
c6_agcrs <- meta %>% filter(tissue == 'Primary' & celltype == 'Fibroblasts') %>% filter(stage %in% c('III', 'IV')) %>% select(c('sampleid', 'AGCRS'))
c6_agcrs$AGCRS <- rescale(c6_agcrs$AGCRS, to = c(0,1))

c6_agcrs_median <- c6_agcrs %>% group_by(sampleid) %>% summarise_each(funs =  median) %>% column_to_rownames('sampleid')

sample_info <- meta %>% select(c('sampleid', 'stage')) %>% unique() %>% rownames_to_column() %>% column_to_rownames('sampleid')
annotation_col <- data.frame(row.names = rownames(c6_agcrs_median), Stage = sample_info[rownames(c6_agcrs_median), 'stage'])
ann_colors <- list(Stage = c('III' = '#4DBBD5FF', 'IV' = '#E64B35FF'))
pdf('Fig_4h.pdf', height = 3, width = 8)
pheatmap(t(c6_agcrs_median), cutree_cols = 2, scale = 'none',
         cluster_rows = F, cluster_cols = T, clustering_method = 'ward.D',
         annotation_colors = ann_colors, annotation_col = annotation_col,
         border_color = 'grey60',
         color = colorRampPalette(colors = rev(brewer.pal(name = 'RdBu', n = 11)))(50),
         cellheight = 20, cellwidth = 10,
         show_colnames = F
)


cluster_info <- dist(c6_agcrs_median) %>% hclust(., method = 'ward.D')
sample_order <- rownames(c6_agcrs_median)[cluster_info$order]

c6_agcrs_median$group <- 'AGCRS_low';
c6_agcrs_median[tail(sample_order, 12), 'group'] <- 'AGCRS_high'
c6_agcrs$group <- factor(c6_agcrs$sampleid, levels = rownames(c6_agcrs_median),
                         labels = c6_agcrs_median$group)

ggplot(c6_agcrs, aes(x = `AGCRS`, y = group, fill = ..density..)) +
  geom_density_ridges_gradient() +
  scale_fill_gradientn(colors = rev(brewer.pal(name = 'Spectral', n = 11))) +
  labs(x = 'AGCRS score', y = '') +
  scale_y_discrete(expand = expansion(add = 0.1, mult = 0.1)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'))
dev.off()


#==========Fig 6b==========
epi_meta <- meta %>% select(c('tissue', 'epitype')) %>% na.omit()
epi_prop <- table(epi_meta$epitype, epi_meta$tissue) %>% prop.table(., 2) %>% as.data.frame()
pdf('fig_4i.pdf', height = 7, width = 6)
ggplot(epi_prop, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = 'identity', position = 'fill', width = 0.8) +
  labs(x = '', y = 'Percentage(%)') +
  scale_y_continuous(expand = c(0.01,0.01), labels = seq(0,100,25)) +
  scale_fill_manual(values = tissue_col) +
  theme_prism(border = T, base_line_size = 0.6, base_fontface = 'bold') +
  theme(legend.position = 'top', legend.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

#==========Supplementary Fig 5d================
pdf('fig_s5d.pdf', height = 5, width = 4.5)
meta %>% filter(celltype %in% c('Epithelials', 'Endocrines', 'Ehdothelials',
                                'Mast cells', 'Macro/Mono', 'T cells', 'B cells', 'Plasmas')) %>%
  ggplot(., aes(x = tissue, y = AGCRS, fill = tissue)) + 
  geom_boxplot(outlier.alpha = 0) +
  facet_wrap(~celltype, nrow = 1) +
  theme_prism(border = T, base_line_size = 0.4, base_size = 12, axis_text_angle = 45) +
  guides(fill = 'none') + labs(x = '', y = 'AGCRS Score') +
  scale_fill_manual(values = tissue_col)
dev.off()

#==========Supplementary Fig 5g==========
pdf('fig_s5g.pdf', height = 7, width = 4)
fibro_prop %>% filter(Var1 %in% c('Peri_C1', 'Peri_C2', 'Peri_C3', 'Mesothelial')) %>%
  ggplot(., aes(x = Var1, y = Freq, fill = Var2,)) +
  geom_bar(stat = 'identity', position = 'fill', width = 0.7, color = 'black') +
  labs(x = '', y = 'Percentage(%)') +
  scale_fill_manual(values = tissue_col) +
  theme_prism(border = T, base_line_size = 0.4, base_fontface = 'bold') +
  scale_y_continuous(expand = c(0.01,0.01), labels = seq(0,100,25)) +
  theme(legend.position = 'top', legend.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

#==========Supplementary Fig 5h==========
dat <- meta %>% filter(celltype %in% c('Pericytes', 'Mesothelials'))
dat$AGCRS <- rescale(dat$AGCRS, to = c(0, 1))
pdf('fig_3i.pdf', height = 5, width = 9)
ggplot(dat, aes(x = tissue, y = AGCRS, fill = tissue)) + 
  geom_boxplot(outlier.alpha = 0, position = position_dodge(0.8), width = 0.5) +
  facet_grid(~fibrotype, scales = 'free_x', space = 'free_x') +
  guides(fill = 'none') + labs(x = '', y = 'AGCRS Score') +
  scale_fill_manual(values = tissue_col) +
  theme_prism(border = T, base_size = 12, axis_text_angle = 45)
dev.off()

