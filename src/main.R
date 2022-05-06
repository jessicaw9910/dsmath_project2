library(magrittr)
library(ggplot2)
library(ggrepel)
library(reshape)
library(dplyr)
library(scales)
library(stringr)
library(svglite)
library(janitor)

# df_drug_info <- read.csv('../assets/drug_data.csv', stringsAsFactors = FALSE)
df_cell_info <- read.csv('../assets/cell_list.csv', stringsAsFactors = FALSE)

# df_drug <- read.csv('../assets/drug.csv', stringsAsFactors = FALSE)
df_cell <- read.csv('../assets/cell.csv', stringsAsFactors = FALSE)

# df_drug_lrm <- read.csv('../assets/drug_lrm.csv', stringsAsFactors = FALSE)
# df_cell_lrm <- read.csv('../assets/cell_lrm.csv', stringsAsFactors = FALSE)

###############################
#### Shorten drug pathways ####
###############################

# df_path <- data.frame(pathway_name = unique(df_drug_info$pathway_name))
# df_path$pathway_short <- c('Pro stab/deg', 'PI3K/MTOR', 'DNA rep',
#                            'EGFR', 'Mitosis', 'Cell cycle', 'Other',
#                            'RTK', 'Apoptosis', 'Unclassified', 'MAPK',
#                            'Misc kinase', 'WNT', 'Genome', 'Hist meth',
#                            'Chromatin other', 'Metabolism', 'p53',
#                            'Cytoskeleton', 'Hormone', 'Hist acet', 'IGF1R',
#                            'ABL', 'JNK/p38')

#########################
#### Clean Drug Data ####
#########################

## create index for later use
df_drug$idx <- seq(1, length(df_drug$drug))

## two entries require manual fixes
idx <- which(df_drug$drug %in% "KRAS\nG12C Inhibitor-12\n1855")
df_drug$drug[idx] <- "KRAS (G12C) Inhibitor-12\n1855"
idx <- which(df_drug$drug %in% "Nutlin-3a\n-\n1047")
df_drug$drug[idx] <- "Nutlin-3a (-)\n1047"

## extract drug name alone
df_drug$drug_name <- sapply(df_drug$drug,
                            function(x) gsub("[\r\n].*", "", x))

## create column with log10(1/IC50)
df_drug$inv_IC50 <- log10(1/exp(df_drug$lnIC50_mean))

## start cluster at 1
knn_col <- colnames(df_drug)[grep(pattern = 'knn', colnames(df_drug))]
for (col in knn_col){
  df_drug[, col] <- df_drug[, col] + 1
}

## merge with drug pathway info
df_drug <- merge(df_drug, df_drug_info[, c('drug_name', 'pathway_name')], 
                 by = 'drug_name', all.x = TRUE)

## remove duplicates
df_drug <- df_drug[-which(duplicated(df_drug$drug)), ]

#########################
#### Clean Cell Data ####
#########################

df_cell_info <- clean_names(df_cell_info)

df_cell_info$tcga_classification[(df_cell_info$name == 'NCCIT')] <- 'UNCLASSIFIED'

df_cell$idx <- seq(1, length(df_cell$cell))

## create column with log10(1/IC50)
df_cell$inv_IC50 <- log10(1/exp(df_cell$lnIC50_mean))

## start cluster at 1
knn_col <- colnames(df_cell)[grep(pattern = 'knn', colnames(df_cell))]
for (col in knn_col){
  df_cell[, col] <- df_cell[, col] + 1
}

df_cell <- merge(df_cell, df_cell_info[, c('name', 'tcga_classification', 'tissue')],
                 by.x = 'cell', by.y = 'name', all.x = TRUE, all.y = FALSE)

df_cell <- df_cell[-which(duplicated(df_cell$cell)), ]

############################
#### Drug K-Means Plots ####
############################

n_all <- length(unique(df_drug$knn_all))
n_pca <- length(unique(df_drug$knn_pca))
n_svd <- length(unique(df_drug$knn_svd))

## plot KNN results for all data
df_drug %>% ggplot(aes(x=tsne1, y=tsne2, 
                       color=pathway_name, shape=factor(knn_all))) +
  geom_point(aes(size = inv_IC50)) +
  geom_text_repel(size = 3, aes(label = toupper(drug_name)), 
                  point.padding = unit(0.3, "lines"), show.legend = FALSE) +
  scale_shape_manual(values = 0:n_all, guide="none") +
  labs(title = "Relative size indicates log10(1/Average IC50)") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_size(guide="none") +
  guides(color=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))
  
## plot KNN results for PCA
df_drug %>% ggplot(aes(x=tsne1, y=tsne2, 
                       color=pathway_name, shape=factor(knn_pca))) +
  geom_point(aes(size = inv_IC50)) +
  geom_text_repel(size = 3, aes(label = toupper(drug_name)), 
                  point.padding = unit(0.3, "lines"), show.legend = FALSE) +
  scale_shape_manual(values = 0:n_pca, guide="none") +
  labs(title = "Relative size indicates log10(1/Average IC50)") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_size(guide="none") +
  guides(color=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))

## plot KNN results for SVD
df_drug %>% ggplot(aes(x=tsne1, y=tsne2, 
                       color=pathway_name, shape=factor(knn_svd))) +
  geom_point(aes(size = inv_IC50)) +
  geom_text_repel(size = 3, aes(label = toupper(drug_name)), 
                  point.padding = unit(0.3, "lines"), show.legend = FALSE) +
  scale_shape_manual(values = 0:n_svd, guide="none") +
  labs(title = "Relative size indicates log10(1/Average IC50)") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  scale_size(guide="none") +
  guides(color=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))

#############################
#### Low-Rank Drug Model ####
#############################
  
colnames(df_drug_lrm) <- paste(rep('V', 10), 
                               seq(1, dim(df_drug_lrm)[2]), sep = "")

df_drug_top10 <- as.data.frame(apply(df_drug_lrm, 2, 
                                     function(x) order(x, decreasing = TRUE)))
df_drug_top10$rank <- seq(1, dim(df_drug_top10)[1])
df_drug_top10 <- melt(df_drug_top10, id.vars = 'rank')
colnames(df_drug_top10) <- c('rank', 'V', 'idx')
df_drug_top10 <- merge(df_drug_top10, df_drug[, c('idx', 'inv_IC50', 'drug_name', 'pathway_name')],
                       by = 'idx', all.x = TRUE)

idx_top <- which(df_drug_top10$rank <= 10)
idx_bottom <- which(df_drug_top10$rank >= 189)

## create dataframe with top and bottom 10 values
df_drug_top <- df_drug_top10[idx_top, ]
df_drug_top$facet <- "Top 10"
df_drug_bottom <- df_drug_top10[idx_bottom, ]
df_drug_bottom$facet <- "Bottom 10"
df_drug_top10_combo <- rbind(df_drug_top, df_drug_bottom)
df_drug_top10_combo$facet <- relevel(as.factor(df_drug_top10_combo$facet), 'Top 10')

ggplot(df_drug_top10_combo, aes(x=rank, y=V, size=inv_IC50, color=pathway_name)) +
  geom_point() +
  facet_wrap(~facet, scales = 'free')

ggplot(df_drug_top10_combo, aes(x=V, y=inv_IC50)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour=pathway_name, size=inv_IC50), alpha=0.5) +
  facet_wrap(~facet, dir = 'v') +
  xlab('Factor') +
  ylab('log10(1/Average IC50)') + 
  scale_size(guide="none") +
  guides(color=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))

ggsave('plots/Top10_Boxplot_Drug.png', width = 12, height = 8)
ggsave('plots/Top10_Boxplot_Drug.svg', width = 12, height = 6)

# df_stacked_top <- df_drug_top10[idx_top, ] %>% 
#   count(V, pathway_name, sort = TRUE)
# df_stacked_bottom <- df_drug_top10[idx_bottom, ] %>% 
#   count(V, pathway_name, sort = TRUE)

df_stacked_top <- df_drug_top10[idx_top, ] %>%
  group_by(V, pathway_name) %>%
  summarise(all_names = paste(drug_name, collapse = ", "))
df_stacked_top$n <- str_count(df_stacked_top$all_names, ",") + 1
df_stacked_top$order <- "Top 10"

df_stacked_bottom <- df_drug_top10[idx_bottom, ] %>%
  group_by(V, pathway_name) %>%
  summarise(all_names = paste(drug_name, collapse = ", "))
df_stacked_bottom$n <- str_count(df_stacked_bottom$all_names, ",") + 1
df_stacked_bottom$order <- "Bottom 10"

# ggplot(df_stacked_top, aes(x = V, y = n, fill = pathway_name)) + 
#   geom_bar(stat = "identity")
# ggplot(df_stacked_bottom, aes(x = V, y = n, fill = pathway_name)) + 
#   geom_bar(stat = "identity")

df_stacked <- rbind(df_stacked_top, df_stacked_bottom)
df_stacked$percent <- df_stacked$n / 10
df_stacked <- df_stacked[order(df_stacked$n), ]
df_stacked <- merge(df_stacked, df_path, by = 'pathway_name', all.x = TRUE)
df_stacked$order <- relevel(as.factor(df_stacked$order), 'Top 10')

sum(df_stacked$percent[which(df_stacked$order == "Top 10")] >= .5)
sum(df_stacked$percent[which(df_stacked$order == "Bottom 10")] >= .5)

ggplot(df_stacked, aes(x = V, y = percent, fill = pathway_name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~order, dir = 'v') +
  geom_text(aes(label=stringr::str_wrap(all_names, 5)), size = 3, position = position_fill(vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('Factor') +
  ylab('Percent') + 
  guides(fill=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))

ggsave('plots/Top10_Stacked Bar_Drug.png', width = 12, height = 8)
ggsave('plots/Top10_Stacked Bar_Drug.svg', width = 12, height = 6)

# ggplot(df_stacked, aes(x = V, y = percent, fill = pathway_name)) + 
#   geom_bar(stat = "identity") +
#   facet_wrap(~order, dir = 'v') +
#   geom_text(aes(label=paste(pathway_short, n, sep = ': ')), size = 3, position = position_fill(vjust = 0.5)) +
#   scale_y_continuous(labels = scales::percent_format()) +
#   xlab('Factor') +
#   ylab('Percent') + 
#   scale_fill_discrete(guide="none")

ggplot(df_stacked, aes(x = V, y = percent, fill = pathway_name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~order, dir = 'v') +
  geom_text(aes(label=paste(pathway_short, percent(percent, accuracy = 1), sep=": ")),
            size = 3, position = position_fill(vjust = 0.5)) +
  theme(legend.position = "none") 

ggsave('../src/plots/Top10_Stacked Bar_MOA.svg', width = 12, height = 6)

df_stacked$all_names[which(df_stacked$n == 10)]
df_stacked$all_names[which(df_stacked$n == 7)]

sum(df_stacked$percent[which(df_stacked$order == "Top 10")] >= .5)
sum(df_stacked$percent[which(df_stacked$order == "Bottom 10")] >= .5)

#############################
#### Low-Rank Cell Model ####
#############################

colnames(df_cell_lrm) <- paste(rep('U', 10),
                               seq(1, dim(df_cell_lrm)[2]), sep = "")

df_cell_top10 <- as.data.frame(apply(df_cell_lrm, 2,
                                     function(x) order(x, decreasing = TRUE)))
df_cell_top10$rank <- seq(1, dim(df_cell_top10)[1])
df_cell_top10 <- melt(df_cell_top10, id.vars = 'rank')
colnames(df_cell_top10) <- c('rank', 'U', 'idx')
df_cell_top10 <- merge(df_cell_top10, df_cell[, c('idx', 'cell', 'inv_IC50', 'tcga_classification', 'tissue')],
                       by = 'idx', all.x = TRUE)

idx_top <- which(df_cell_top10$rank <= 10)
idx_bottom <- which(df_cell_top10$rank >= 800)

## create dataframe with top and bottom 10 values
df_cell_top <- df_cell_top10[idx_top, ]
df_cell_top$facet <- "Top 10"
df_cell_bottom <- df_cell_top10[idx_bottom, ]
df_cell_bottom$facet <- "Bottom 10"
df_cell_top10_combo <- rbind(df_cell_top, df_cell_bottom)
df_cell_top10_combo$facet <- relevel(as.factor(df_cell_top10_combo$facet), 'Top 10')

ggplot(df_cell_top10_combo, aes(x=U, y=inv_IC50)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour=tissue, size=inv_IC50), alpha=0.5) +
  facet_wrap(~facet, dir = 'v') +
  xlab('Factor') +
  ylab('log10(1/Average IC50)') + 
  scale_size(guide="none") +
  guides(color=guide_legend(title="Tissue", ncol=1, 
                            override.aes = list(size=6)))

ggsave('plots/Top10_Boxplot_Cell.png', width = 12, height = 8)
ggsave('plots/Top10_Boxplot_Cell.svg', width = 12, height = 6)

df_stacked_top <- df_cell_top10[idx_top, ] %>%
  group_by(U, tissue) %>%
  summarise(all_names = paste(cell, collapse = ", "))
df_stacked_top$n <- str_count(df_stacked_top$all_names, ",") + 1
df_stacked_top$order <- "Top 10"

df_stacked_bottom <- df_cell_top10[idx_bottom, ] %>%
  group_by(U, tissue) %>%
  summarise(all_names = paste(cell, collapse = ", "))
df_stacked_bottom$n <- str_count(df_stacked_bottom$all_names, ",") + 1
df_stacked_bottom$order <- "Bottom 10"

df_stacked <- rbind(df_stacked_top, df_stacked_bottom)
df_stacked$percent <- df_stacked$n / 10
df_stacked <- df_stacked[order(df_stacked$n), ]
df_stacked$order <- relevel(as.factor(df_stacked$order), 'Top 10')

sum(df_stacked$percent[which(df_stacked$order == "Top 10")] >= .5)
sum(df_stacked$percent[which(df_stacked$order == "Bottom 10")] >= .5)

ggplot(df_stacked, aes(x = U, y = percent, fill = tissue)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~order, dir = 'v') +
  geom_text(aes(label=stringr::str_wrap(all_names, 5)), size = 2.5, position = position_fill(vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('Factor') +
  ylab('Percent') + 
  guides(fill=guide_legend(title="Tissue", ncol=1, 
                           override.aes = list(size=6)))

ggsave('plots/Top10_Stacked Bar_Cell.png', width = 12, height = 8)
ggsave('plots/Top10_Stacked Bar_Cell.svg', width = 12, height = 6)

#################################
#### Drug K-Means Clustering ####
#################################

## for drug clusters
cols <- c('drug_name', 'knn_all', 'knn_pca', 'knn_svd', 'pathway_name')

df_knn_drug <- df_drug[, cols]
df_knn_drug <- melt(df_knn_drug, id.vars = c('drug_name', 'pathway_name'))
colnames(df_knn_drug)[dim(df_knn_drug)[2]] <- 'cluster'
colnames(df_knn_drug)[dim(df_knn_drug)[2]-1] <- 'representation'

df_knn_drug <- df_knn_drug %>%
  group_by(cluster, pathway_name, representation) %>%
  summarise(all_names = paste(drug_name, collapse = ", "))

df_knn_drug$count <- str_count(df_knn_drug$all_names, ",") + 1
df_knn_drug$cluster <- as.factor(df_knn_drug$cluster)

ggplot(df_knn_drug, aes(x=cluster, y=count, fill=pathway_name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v')

df_sum_drug <- aggregate(count ~ cluster + representation, df_knn_drug, sum)
colnames(df_sum_drug)[dim(df_sum_drug)[2]] <- 'total'

df_knn_drug <- merge(df_knn_drug, df_sum_drug, by=c('cluster', 'representation'))

df_knn_drug$percent <- df_knn_drug$count / df_knn_drug$total
df_knn_drug <- merge(df_knn_drug, df_path, by='pathway_name')
df_knn_drug$representation <- gsub(pattern = 'knn_', replacement = '', df_knn_drug$representation)
df_knn_drug$representation <- sapply(df_knn_drug$representation, function(x) ifelse(x!='all', toupper(x), 'All'))

# df_knn_drug$dominant <- (df_knn_drug$percent >= 0.4)
# df_knn_drug[df_knn_drug$dominant == TRUE, ] %>% View()
# table(df_knn_drug$representation[df_knn_drug$dominant == TRUE])
# df_knn_drug$label <- mapply(FUN = function(x, y, z) ifelse(x == TRUE, y , z),
#                             df_knn_drug$dominant, 
#                             df_knn_drug$pathway_name, 
#                             df_knn_drug$cluster)

df_temp <- df_knn_drug[(df_knn_drug$percent >= 0.4), c('pathway_name', 'cluster', 'representation')]
colnames(df_temp)[1] <- 'label'

df_knn_drug <- merge(df_knn_drug, df_temp, by = c('cluster', 'representation'), all.x = TRUE)

ggplot(df_knn_drug[!is.na(df_knn_drug$label), ], aes(x=label, y=percent, fill=pathway_short)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v') +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label=stringr::str_wrap(all_names, 5)),
            size = 2, position = position_fill(vjust = 0.5)) +
  # geom_text(aes(label=paste(pathway_short, count, sep=": ")),
  #           size = 3, position = position_fill(vjust = 0.5)) +
  xlab('Cluster') +
  ylab('Percent') +
  # theme(legend.position = "none")
  guides(fill=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))

ggsave('plots/Drug_Clusters.png', width = 12, height = 8)
ggsave('plots/Drug_Clusters.svg', width = 12, height = 6)

## unlabeled
df_unlabeled_drug <- df_knn_drug[is.na(df_knn_drug$label), ]

ggplot(df_unlabeled_drug, aes(x=cluster, y=percent, fill=pathway_short)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v') +
  scale_y_continuous(labels = scales::percent_format()) +
  # geom_text(aes(label=stringr::str_wrap(all_names, 5)),
  #           size = 2, position = position_fill(vjust = 0.5)) +
  geom_text(aes(label=paste(pathway_short, count, sep=": ")),
            size = 2, position = position_fill(vjust = 0.5)) +
  xlab('Cluster') +
  ylab('Percent') +
  theme(legend.position = "none")
  # guides(color=guide_legend(title="Drug Pathway", ncol=1, 
  #                           override.aes = list(size=6)))

#################################
#### Cell K-Means Clustering ####
#################################

create_clustdf <- function(df, cols_total, cols_melt, cols_group, col_sum){
  df_output <- df[, cols_total]
  df_output <- melt(df_output, id.vars = cols_melt)
  colnames(df_output)[dim(df_output)[2]] <- 'cluster'
  colnames(df_output)[dim(df_output)[2]-1] <- 'representation'

  df_output <- df_output %>%
    group_by_at(cols_group) %>%
    summarise(all_names = paste(!! sym(col_sum), collapse = ", "))
  
  df_output$count <- str_count(df_output$all_names, ",") + 1
  df_output$cluster <- as.factor(df_output$cluster)
  
  df_sum <- aggregate(count ~ cluster + representation, df_output, sum)
  colnames(df_sum)[dim(df_sum)[2]] <- 'total'
  
  df_output <- merge(df_output, df_sum, by=c('cluster', 'representation'))
  
  df_output$percent <- df_output$count / df_output$total
  df_output$representation <- gsub(pattern = 'knn_', replacement = '', df_output$representation)
  df_output$representation <- sapply(df_output$representation, function(x) ifelse(x!='all', toupper(x), 'All'))
  
  return(df_output)
}


cols <- c('cell', 'knn_all', 'knn_knn', 'knn_elow', 
          'tissue', 'tcga_classification')

df_knn_tissue <- create_clustdf(df_cell,
                                cols_total = c('cell', 'knn_all', 'knn_knn', 'knn_elow', 'tissue'),
                                cols_melt = c('cell', 'tissue'),
                                cols_group = c('cluster', 'tissue', 'representation'),
                                col_sum = 'cell')

df_knn_tcga <- create_clustdf(df_cell,
                              cols_total = c('cell', 'knn_all', 'knn_pca', 'knn_svd', 'tcga_classification'),
                              cols_melt = c('cell', 'tcga_classification'),
                              cols_group = c('cluster', 'tcga_classification', 'representation'),
                              col_sum = 'cell')

ggplot(df_knn_cell, aes(x=cluster, y=count, fill=tcga_classification)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v')

ggplot(df_knn_tcga, aes(x=cluster, y=count, fill=tcga_classification)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v')

df_temp <- df_knn_tissue[(df_knn_tissue$percent >= 0.4), c('tissue', 'cluster', 'representation')]
colnames(df_temp)[1] <- 'label'
df_knn_tissue <- merge(df_knn_tissue, df_temp, by = c('cluster', 'representation'), all.x = TRUE)

ggplot(df_knn_tissue[!is.na(df_knn_tissue$label), ], aes(x=label, y=percent, fill=tissue)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label=stringr::str_wrap(paste(all_names, count, sep = ' - '), width = 53)),
            size = 2.5, position = position_fill(vjust = 0.5)) +
  # geom_text(aes(label=paste(pathway_short, count, sep=": ")),
  #           size = 3, position = position_fill(vjust = 0.5)) +
  xlab('Cluster') +
  ylab('Percent') +
  # theme(legend.position = "none")
  guides(fill=guide_legend(title="Tissue", ncol=1, 
                           override.aes = list(size=6)))

ggsave('plots/Tissue_Clusters.png', width = 12, height = 8)
ggsave('plots/Tissue_Clusters.svg', width = 12, height = 6)

## unlabeled
df_unlabeled_tissue <- df_knn_cell[is.na(df_knn_tissue$label), ]

ggplot(df_unlabeled_cell, aes(x=cluster, y=percent, fill=tissue)) +
  geom_bar(stat = "identity") +
  facet_wrap(~representation, dir = 'v') +
  scale_y_continuous(labels = scales::percent_format()) +
  # geom_text(aes(label=stringr::str_wrap(all_names, 5)),
  #           size = 2, position = position_fill(vjust = 0.5)) +
  geom_text(aes(label=paste(tissue, count, sep=": ")),
            size = 2, position = position_fill(vjust = 0.5)) +
  xlab('Cluster') +
  ylab('Percent') +
  theme(legend.position = "none")
# guides(color=guide_legend(title="Drug Pathway", ncol=1, 
#                           override.aes = list(size=6)))
