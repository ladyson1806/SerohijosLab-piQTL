library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(Rtsne)


#### Config
base_dir <- getwd()

#### Load PPI matrix
ppi_matrix <- "/../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_noPPI-labelled_for_rMVP.csv"  # ToDo: modify input path
ppi_file_name <- paste0(base_dir, ppi_matrix)
ppi_master <- read_csv(ppi_file_name)

sample_table <- "/../results/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_noPPI-labelled_metadata.csv"  # ToDo: modify input path
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
tsne_plot_ppi <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$PPI)
t1 <- ggplot(tsne_plot_ppi) + geom_point(aes(x=x, y=y, color=col)) + guides(colour= "none")
tsne_plot_drug <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$Drug)
t2 <- ggplot(tsne_plot_drug) + geom_point(aes(x=x, y=y, color=col)) + guides(colour= "none")
tsne_plot_mtx <- data.frame(x = ppi_tsne_out$Y[,1], y = ppi_tsne_out$Y[,2], col = ppi_tsne$MTX)
t3 <- ggplot(tsne_plot_mtx) + geom_point(aes(x=x, y=y, color=col)) + guides(colour= "none")

setEPS()
postscript(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_noPPI-labelled_logratio_Fitness_tnsne.eps", width = 12, height = 6, family = "Helvetica", pointsize=12)  # ToDo: modify output path
# print the result (in this case a ggplot)
grid.arrange(t1, t2, t3, nrow =1)
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
postscript(file = "../manuscript/figures/EXT_FIGURE_4//minus_ref_noPPI-labelled_logratio_Fitness_PCs_per_PPI.eps", width = 12, height = 6, family = "Helvetica", pointsize=12)  # ToDo: modify output path
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
  geom_point()

p2_DRUG <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = Drug)) +
  geom_point()


p3_DRUG <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC2, y = PC3, colour = Drug)) +
  geom_point()

setEPS()
postscript(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_noPPI-labelled_logratio_Fitness_PCs_per_Drug.eps", width = 18, height = 6, family = "Helvetica", pointsize=12)  # ToDo: modify output path
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
  geom_point()

p2_MTX <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = MTX)) +
  geom_point()


p3_MTX <- p$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "label") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "label") %>%
  # create the plot
  ggplot(aes(x = PC2, y = PC3, colour = MTX)) +
  geom_point()

setEPS()
postscript(file = "../manuscript/figures/EXT_FIGURE_4/minus_ref_noPPI-labelled_logratio_Fitness_PCs_per_MTX.eps", width = 18, height = 6, family = "Helvetica", pointsize=12)  # ToDo: modify output path
grid.arrange(p1_MTX, p2_MTX, p3_MTX, nrow = 1)
dev.off()
