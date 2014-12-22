# Simple visualization of metagenomic sample principal components
# computed by RTG similarity command.

# Requires: R (version 2.15 or later), rgl (will attempt to install)
#
# Version: 1.0


# The similarity.pca file in this package has been output by the
# similarity command. The first three columns are the principal
# component values, and the fourth is the name of the sample. In this
# example, if we consider the sample name to be composed of four
# sub-fields, the third of these indicates the site from which the
# sample was taken. In this visualization, we use the site to assign a
# color. To run the visualization, start R from this directory and
# issue the command:
#
# source("similarity-simple.R")


# Attempt to install the rgl package if it is not already available.
# Once the package has been installed it will be available for future
# runs, so this cost is only incurred the first time the script is used.
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if (!is.installed("rgl")) {
  cat("Attempting to install rgl package.  This may take a while or fail.\n")
  cat("For extra hints look at the comments in the script.\n")
  install.packages("rgl")
}

# Required library
library(rgl)

# Load the principal component analysis file created by RTG similarity
pca <- read.table("similarity.pca", FALSE, "\t", "")
colnames(pca) <- c("X", "Y", "Z", "Sample")

# Make a color lookup table to be used with the sample site
colorlookup <- c('AN', 'BM', 'PF', 'S', 'SP', 'TD')

# Split sample names by _, and use the third field to assign a color id from the lookup table
color <- match(sapply(strsplit(as.character(pca$Sample), "_"), "[[", 3), colorlookup)
pca <- data.frame(pca, color);


# Plot the samples in PCA space, using the assigned colors
open3d(windowRect=c(100,100,1024,768))
plot3d(pca$X, pca$Y, pca$Z, xlab="A", ylab="B", zlab="C", col=pca$color, size=1, type='s')
text3d(pca$X, pca$Y, pca$Z, pca$Sample, adj=-0.15, col=pca$color, cex=0.75)
