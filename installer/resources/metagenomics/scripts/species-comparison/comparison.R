#!/usr/bin/Rscript

# Example script for the comparison of outputs from the RTG species module.
# It takes two species output directories and produces a visual comparison
# of the two samples.
#
# Requires: R (version 2.14 or later), ggplot2 (will attempt to install)
# Author: Mehul Rathod <mehul@realtimegenomics.com>
# Author: Sean A. Irvine <sean@realtimegenomics.com>
# Version: 1.0

# Parse the command line arguments.  If they are not what is expected then
# produce the usage message.
args <- commandArgs(TRUE)
if (length(args) != 4) {
  cat("USAGE: comparison.R species-output-dir1 label1 species-output-dir2 label2\n")
  cat("\tspecies-output-dir1   first sample to compare\n")
  cat("\tlabel1                arbitrary label for first sample\n")
  cat("\tspecies-output-dir2   second sample to compare\n")
  cat("\tlabel2                arbitrary label for second sample\n")
  q(status=1)
}

# For best performance you should choose a mirror of CRAN close to you.
# See http://cran.r-project.org/mirrors.html for complete list of mirrors.
options("repos"="http://cran.us.r-project.org")    ## US mirror
#options("repos"="http://cran.csiro.au/ ") ## Australian mirror
#options("repos"="http://cran.stat.auckland.ac.nz/") ## NZ mirror
#options("repos"="http://www.stats.bris.ac.uk/R/") ## UK mirror

# Attempt to install the ggplot2 package if it is not already available.
# Once the package has been installed it will be available for future
# runs, so this cost is only incurred the first time the script is used.
#
# The package can be installed manually from the command line with
#   wget http://cran.r-project.org/src/contrib/ggplot2_0.9.3.tar.gz
#   R CMD INSTALL ggplot2_0.9.3.tar.gz
# but in this approach you will have to sort out many dependencies.
#
# Installation can fail for a variety of reasons:
#  - ggplot2 requires R 2.14 or greater
#  - permission problem writing to library directory
#    (can be circumvented by setting R_LIBS environment variable)
#  - lack of Internet connectivity

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if (!is.installed("ggplot2")) {
  cat("Attempting to install ggplot2 package.  This may take a while or fail.\n")
  cat("For extra hints look at the comments in the script.\n")
  install.packages("ggplot2")
}

library(ggplot2)

# This completes the set up.

# We are now ready to retrieve the species information from the directories
# provided on the command line.  We assume these files are exactly as
# produced by the RTG species module.  In particular, we look for a file
# called "species.tsv" which is assumed to be sorted by the abundance.
# Although "species.tsv" does have a header, we skip it here because it
# look just like all the other comment lines at the top.

file1name <- paste(args[1], "/species.tsv", sep="")
sample1 <- read.csv(file1name, comment.char="#", header=F, sep="\t")
label1 <- args[2]

file2name <- paste(args[3], "/species.tsv", sep="")
sample2 <- read.csv(file2name, comment.char="#", header=F, sep="\t")
label2 <- args[4]

# Attach column names to the data frames. Use "_" rather than "-" in the
# column names so that R column references work.
header <- c("abundance",
            "abundance_low",
            "abundance_high",
            "DNA_fraction",
            "DNA_fraction_low",
            "DNA_fraction_high",
            "confidence",
            "coverage_depth",
            "coverage_breadth",
            "mapped_reads",
            "has_reference",
            "taxa_count",
            "taxon_id",
            "parent_id",
            "rank",
            "taxonomy_name")
colnames(sample1) <- header
colnames(sample2) <- header

# We only want a comparison between "species" have sequences associated with
# them.  This prevents the comparison from getting cluttered from higher
# taxonomic items.  This is controlled by the has_sequence column in the
# "species.tsv" file.
sample1 <- sample1[sample1$has_reference == "Y",]
sample2 <- sample2[sample2$has_reference == "Y",]

# Names of the columns we need to retain to make the plot and determine
# the ordering of the species.
retained <- c("taxonomy_name", "abundance", "abundance_low", "abundance_high")

# Merge the samples replacing any missing values with 0.  Compute the max
# of the abundances from the two samples, finally limit to the top 20 rows.
merge1 <- sample1[retained]
merge2 <- sample2[retained]
merged <- merge(merge1, merge2, by="taxonomy_name", all.x=T, all.y=T)
merged[is.na(merged)] <- 0
merged$max_abundance <- with(merged, pmax(abundance.x, abundance.y))
merged <- merged[order(-merged$max_abundance),]
merged <- merged[1:min(20, nrow(merged)),]

# Resplit into the original samples based on those retained
sample1 <- merged[c("taxonomy_name", "abundance.x", "abundance_low.x", "abundance_high.x", "max_abundance")]
sample2 <- merged[c("taxonomy_name", "abundance.y", "abundance_low.y", "abundance_high.y", "max_abundance")]
# Reset column names in both samples
revised_header <- c("taxonomy_name", "abundance", "abundance_low", "abundance_high", "max_abundance")
colnames(sample1) <- revised_header
colnames(sample2) <- revised_header

# Add the label to each sample, based on the command line argument
sample1$Sample <- label1
sample2$Sample <- label2

# Combine the two samples vertically into a single frame in preparation for plotting
combined <- rbind(sample1, sample2)

# Uncomment the following lines, if you want the table of plotted data to be
# printed on the terminal before graphing.
#options(width=10000)
#combined

# Sort the combined data into reverse abundance order.
# Plot the comparison using a log scale for the abundance and with error bars.
plot <- ggplot(data=combined, aes(x=reorder(taxonomy_name, max_abundance), y=abundance, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=abundance_low, ymax=abundance_high), size=.3,  width=.7, position=position_dodge(.9)) +
  xlab("species") +
  ylab("abundance") +
  theme_bw() +
  # Next line circumvents ordering of lengend bug in R when using coord_flip()
  guides(fill = guide_legend(reverse = T)) +
  ggtitle("Comparison of top 20 species") +
  coord_flip()

# Select the name and dimensions of output file. Write the plot to the
# file and then turn of the output device.
png(filename="comparison.png", width = 1024, height = 768)
print(plot)

# The assignment below prevents the message "null device" from
# appearing in the terminal window.
ignore <- dev.off()

# Let the user know where the output is
cat("Wrote graph in \"comparison.png\"\n")
