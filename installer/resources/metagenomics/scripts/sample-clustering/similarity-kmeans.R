# Perform k-means clustering of metagenomic samples from principal
# component analysis outputs computed by the RTG similarity command.

# Requires: R (version 2.15 or later), rgl (will attempt to install),
# fpc (will attempt to install)

# The similarity.pca file in this package has been output by the
# similarity command. The first three columns are the principal
# component values, and the fourth is the name of the sample. In this
# example, we perform k-means clustering within the PCA dimensions,
# and set the color based on the cluster to which each sample is
# assigned.  To run the visualization, start R from this directory and
# issue the command:
#
# source("similarity-kmeans.R")


# Attempt to install the rgl and fpc packages if it is not already available.
# Once the package has been installed it will be available for future
# runs, so this cost is only incurred the first time the script is used.
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if (!is.installed("rgl")) {
  cat("Attempting to install rgl package.  This may take a while or fail.\n")
  cat("For extra hints look at the comments in the script.\n")
  install.packages("rgl")
}
if (!is.installed("fpc")) {
  cat("Attempting to install fpc package.  This may take a while or fail.\n")
  cat("For extra hints look at the comments in the script.\n")
  install.packages("fpc")
}

# Required library
library(fpc)
library(rgl)

# Load the principal component analysis file created by RTG similarity
pca <- read.table("similarity.pca", FALSE, "\t", "")
colnames(pca) <- c("X", "Y", "Z", "Sample")

# Create matrix from PCA values for clustering
x <- c(pca$X, pca$Y, pca$Z)
y <- matrix(x, length(x)/3, 3)

# Compute k-means number of clusters
result <- pamk(y)
fit <- kmeans(y, result$nc, nstart=10)

# Add cluster ids and sample names to data 
y <- data.frame(y, fit$cluster)
y <- data.frame(y, pca$Sample)

# Plot the samples in PCA space, using cluster id for color
open3d(windowRect=c(100,100,1024,768))
plot3d(y$X1, y$X2, y$X3, xlab="A", ylab="B", zlab="C", col=y$fit.cluster, size=1, type='s')
text3d(y$X1, y$X2, y$X3, y$pca.Sample, adj=-0.15, col=y$fit.cluster, cex=0.75)
