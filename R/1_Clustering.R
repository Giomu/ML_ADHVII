################################################################################
########### First Analysis: Dimensionality Reduction and Clustering ############
################################################################################
# This script performs dimensionality reduction and clustering on real data. 
# The dataset contains 116 observations and 13 variables, the same as in the
# Manuscript.
# Comparison is between two methods of dimensionality reduction: UMAP and t-SNE.
# After reducing the dimensionality, a Gaussian Mixture Models (GMM) is applied
# for clustering and the evaluation of clustering results is performed using:
#  - Within-cluster sum of squares (WSS) 
#  - Average silhouette width
# 
# The script produces visualizations for both UMAP and t-SNE projections, 
# as well as clustering results and metrics.
################################################################################

set.seed(983)
# Load the necessary libraries
library(umap)
library(mclust)
library(ggplot2)
library(Rtsne)
library(fpc)

# Load data
load("data/data1.rda")
#load("data/synth1.rda")
head(data1)
rownames(data1) <- data1$ID
data1$ID <- NULL
# Transform the data: log2 transformation and scaling
data1[, c(2:9)] <- log2(data1[, c(2:9)] + 1)
# scale the data
data1[, c(2:13)] <- scale(data1[, c(2:13)], center = T, scale = T)

##############################
# Dimensionality Reduction: UMAP
##############################
set.seed(1778) 
# Apply UMAP on the dataset, excluding the first column (assumed to be Infectious status)
um <- umap::umap(data1[,-c(1)], method = 'umap-learn', preserve.seed = T, 
                 min_dist = 0.5)

# Convert UMAP results into a data frame
umap_df <- data.frame(um$layout)
umap_df$Group <- as.factor(data1$Infection_0_1pre_2post)
umap_df$ID <- rownames(umap_df) 

# Plot UMAP projection
ggplot(umap_df, aes(x = X1, y = X2, color = Group)) + 
  geom_point(size = 3) +
  scale_color_manual(labels = c("No Infection", "Infection pre-boost", 
                                "Infection post-boost"),
                     values = c("#66a182", "#2e4057", "#edae49")) +
  labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Dimensionality Reduction (UMAP)') +
  theme_minimal()

# Clustering with Gaussian Mixture Model (GMM) on UMAP-reduced data (unsupervised)
um_gmm = mclust::Mclust(umap_df[, c(1, 2)]) 
# Print summary of clustering results
summary(um_gmm)
# Plot the clustering results
plot(um_gmm, "classification") 
plot(um_gmm, "density")

##############################
# Dimensionality Reduction: t-SNE
##############################
set.seed(1848)
# Apply t-SNE on the dataset, excluding the first column (assumed to be Infectious status)
tsne <- Rtsne::Rtsne(data1[, -c(1)], perplexity = 37, normalize=FALSE)

# Convert t-SNE results into a data frame
tsne_df <- data.frame(tsne$Y)
tsne_df$Group <- as.factor(data1$Infection_0_1pre_2post)
tsne_df$ID <- rownames(data1)

# Plot t-SNE projection
ggplot(tsne_df, aes(x = X1, y = X2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(labels = c("No Infection", "Infection pre-boost", 
                                "Infection post-boost"),
                     values = c("#66a182", "#2e4057", "#edae49")) +
  labs(x = 't-SNE 1', y = 't-SNE 2', title = 'Dimensionality Reduction (t-SNE)') +
  theme_minimal()

# Clustering with Gaussian Mixture Model (GMM) on t-SNE-reduced data (unsupervised)
t_gmm <- mclust::Mclust(tsne_df[, c(1, 2)])

# Print summary of clustering results
summary(t_gmm)
# Plot the clustering results
plot(t_gmm, "classification")
plot(t_gmm, "density")
# Add clustering results to the t-SNE data frame
tsne_df$Cluster <- as.factor(t_gmm$classification)

##############################
# Clustering Performance Evaluation
##############################
# Compute clustering statistics for UMAP-based clustering
cs_um_gmm <- fpc::cluster.stats(dist(umap_df[1:2]), um_gmm$classification)
stats_um_gmm <- cs_um_gmm[c("within.cluster.ss","avg.silwidth")]

# Compute clustering statistics for t-SNE-based clustering
cs_ts_gmm <- fpc::cluster.stats(dist(tsne_df[1:2]), t_gmm$classification)
stats_t_gmm <- cs_ts_gmm[c("within.cluster.ss","avg.silwidth")]

# Combine statistics and print Comparison
stats <- rbind(stats_um_gmm, stats_t_gmm)
rownames(stats) <- c("UMAP", "t-SNE")
stats <- as.data.frame(stats)
stats

##############################
# Final Clustering Visualization 
##############################
# Plot t-SNE projection with clustering results as the best-performing method
ggplot(tsne_df, aes(x = X1, y = X2, color = Cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#2e4057", "firebrick"),
                     labels = c("Cluster 1", "Cluster 2"),
                     name = "GMM Clusters"
                     ) +
  labs(x = 't-SNE 1', y = 't-SNE 2', title = '') +
  theme_minimal()
