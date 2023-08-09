library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(Rtsne)


#### Config
base_dir <- getwd()

#### Load PPI matrix
ppi_matrix <- "/../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_for_rMVP.csv"
ppi_file_name <- paste0(base_dir, ppi_matrix)
ppi_master <- read_csv(ppi_file_name)

sample_table <- "/../results/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_metadata.csv"
metadata_file_name <- paste0(base_dir, sample_table)
sample_info <- read_csv(metadata_file_name)

ppi_master <- subset(ppi_master, select = -1)
ppi_master_heatmap <- ppi_master %>%
    as.matrix() %>%
    t()
heatmap(ppi_master_heatmap, scale = "column")


ppi_master_pca <- ppi_master %>%
  as.matrix() %>%
  t()
p <- prcomp(ppi_master_pca)

pc_eigenvalues <- p$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)),
                         variance = pc_eigenvalues) %>%
  # add a new column with the percent variance
  mutate(pct = variance / sum(variance) * 100) %>%
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

p0 <- pc_eigenvalues %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

p0

pc_scores <- p$x

pc_scores <- pc_scores %>%
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label")


ppi_tsne  <- pc_scores %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label")

ppi_tsne_out <- Rtsne(pc_scores[2:6])
# tsne_plot_ppi <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$PPI)
# t1 <- ggplot(tsne_plot_ppi) + geom_point(aes(x=x, y=y, color=col)) + guides(colour= "none") + theme_classic()
tsne_plot_drug <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$Drug) 
t2 <- ggplot(tsne_plot_drug) + geom_point(aes(x=x, y=y, color=col)) + labs(x='Dimension 1', y='Dimension 2') + guides(color=guide_legend(ncol=1)) + theme_classic() + theme(legend.text=element_text(size=7), legend.position="bottom", axis.title.x=element_text(size=7), axis.title.y=element_text(size=7), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7)) + coord_fixed()
tsne_plot_mtx <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$MTX)
t3 <- ggplot(tsne_plot_mtx) + geom_point(aes(x=x, y=y, color=col)) + labs(x='Dimension 1', y='Dimension 2') + guides(color=guide_legend(nrow=2)) + theme_classic() + theme(legend.text=element_text(size=7), legend.position="bottom", axis.title.x=element_text(size=7), axis.title.y=element_text(size=7), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7)) + coord_fixed()

setEPS()
postscript(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_logratio_Fitness_tnsne.eps", width = 7.1, height = 3.5, family = "Helvetica", pointsize=7)
# print the result (in this case a ggplot)
grid.arrange(t2, t3, nrow =1)
dev.off()

p1_PPI <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = PPI)) +
  geom_point() +
  guides(colour= "none")

p2_PPI <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = PPI)) +
  geom_point() +
  guides(colour= "none")

p3_PPI <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC2, y = PC3, colour = PPI)) +
  geom_point() +
  guides(colour= "none")

setEPS()
postscript(file = "../manuscript/figures/EXT_FIGURE_4//minus_ref_logratio_Fitness_PCs_per_PPI.eps", width = 12, height = 6, family = "Helvetica", pointsize=12)
# print the result (in this case a ggplot)
grid.arrange(p1_PPI, p2_PPI, p3_PPI,  nrow =1)
dev.off()


p1_DRUG <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = Drug)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
    x = paste0("PC1 ", "(",round(pc_eigenvalues$pct[1], 2),"%)"),
    y = paste0("PC2 ", "(",round(pc_eigenvalues$pct[2], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)
  )

p2_DRUG <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = Drug)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
    x = paste0("PC1 ", "(",round(pc_eigenvalues$pct[1], 2),"%)"),
    y = paste0("PC3 ", "(",round(pc_eigenvalues$pct[3], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)
  )


p3_DRUG <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC2, y = PC3, colour = Drug)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
    x = paste0("PC2 ", "(",round(pc_eigenvalues$pct[2], 2),"%)"),
    y = paste0("PC3 ", "(",round(pc_eigenvalues$pct[3], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)
  )

pdf(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_logratio_Fitness_PCs_per_Drug.pdf", width = 6, height = 2, family = "Helvetica", pointsize=5)
# print the result (in this case a ggplot)
grid.arrange(p1_DRUG, p2_DRUG, p3_DRUG, nrow = 1)
dev.off()


p1_MTX <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = MTX)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
    x = paste0("PC1 ", "(",round(pc_eigenvalues$pct[1], 2),"%)"),
    y = paste0("PC2 ", "(",round(pc_eigenvalues$pct[2], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)
  )

p2_MTX <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = MTX)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
    x = paste0("PC1 ", "(",round(pc_eigenvalues$pct[1], 2),"%)"),
    y = paste0("PC3 ", "(",round(pc_eigenvalues$pct[3], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7)
  )


p3_MTX <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC2, y = PC3, colour = MTX)) +
  geom_point(size=0.2, alpha=0.4) +
  labs(
      x = paste0("PC2 ", "(",round(pc_eigenvalues$pct[2], 2),"%)"),
      y = paste0("PC3 ", "(",round(pc_eigenvalues$pct[3], 2),"%)")
  ) +
  guides(colour= "none") +
  theme_classic() +
  theme(axis.title.x=element_text(size=7),
        axis.title.y=element_text(size=7),
                axis.text.x=element_text(size=7),
                axis.text.y=element_text(size=7)
  )


pdf(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_logratio_Fitness_PCs_per_MTX.pdf", width = 6, height = 2, family = "Helvetica", pointsize=5)
grid.arrange(p1_MTX, p2_MTX, p3_MTX, nrow = 1)
dev.off()
