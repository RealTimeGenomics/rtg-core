# RTG Core Non-commercial

Copyright (c) 2018 Real Time Genomics Ltd

This software is provided under the Real Time Genomics Ltd Software
Licence Agreement for Academic Non-commercial Research Purposes
only. See [LICENSE.txt](LICENSE.txt)

If you require a license for commercial use or wish to purchase
commercial support contact us via info@realtimegenomics.com.

## Introduction

There are many software packages available for next-gen sequencing
analysis, here are some of the reasons to use RTG:

* Comprehensive, cohesive, fast, sensitive and accurate variant
  calling pipeline from FASTQ through to high-quality variants.
* Robust, easy to use, and well supported.
* Fast and accurate alignment, supporting popular sequencing
  technologies, including alignment of Complete Genomics Inc reads.
* Automatic base-quality recalibration computed on-the-fly during
  mapping and applied during variant calling.
* Variant calling both simple SNPs and complex haplotypes,
  automatically realigning when necessary.
* Pedigree-aware pipeline, joint Bayesian calling of multiple samples
  in pedigrees with varying degrees of relatedness, from unrelated,
  trios, quads, half-siblings, larger families, multi-generation pedigrees.
* Pedigree-aware joint calling automatically includes de novo
  discovery, including specific confidence scoring that the variant is
  de novo. Again, this includes information from larger families where
  available, not just trios. Offspring variant calls are automatically
  phased from inheritance.
* Variant calling allows specification of population priors in either
  from general population variants or those in previously called
  samples.
* Sex-aware pipeline, automatically mapping to the correct chromosomes
  for the sex of the sample, performing diploid or haploid variant
  calling as appropriate (including PAR regions).
* Joint Bayesian somatic variant detection, including support for
  contamination, and providing post-calling updated contamination
  estimate. Supports site-specific somatic mutation priors that allow
  databases such as dbSNP/COSMIC/etc. to support the somatic calling.
* Tools for quality control and sample checking: coverage analysis,
  detection of chromosome abnormalities, verification of sex, sample
  mislabelling, incorrect pedigree.
* Includes sophisticated variant comparison tools for benchmarking and
  ROC analysis to guide variant scoring and filtering.  Including tool
  for visualization of ROC information between runs.
* Fast and comprehensive metagenomic analysis pipelines.
* Species frequency composition and abundance analysis, including
  bounds estimates and confidence scores both at the species and
  higher taxon levels.
* Species reported using both tabular data reports and interactive
  Krona HTML pie charts.
* Contaminant filtering and sample similarity analysis.
* Includes tools for building and managing species reference databases
  that include taxonomic structure (prebuilt databases are also
  available).
* Search reads directly against protein databases for functional
  analysis.
* Reproducible results: unlike many other tools, the results don't
  change when you enable multi-threading or when you repeat the same
  run twice.


RTG Core is available pre-packaged directly from our
[website](http://realtimegenomics.com/products/rtg-core/), or follow
the instructions below to build from this repository.

## Support

A user manual is included within the installation in both PDF and HTML
versions. These may also be viewed
online ([HTML](https://realtimegenomics.github.io/rtg-core/index.html),
[PDF](https://cdn.rawgit.com/RealTimeGenomics/rtg-core/master/installer/resources/core/RTGOperationsManual.pdf)).

You can use the commands in RTG Core to format your own reference
datasets, or download common
[pre-formatted references](http://realtimegenomics.com/news/pre-formatted-reference-datasets/)
from our website.

An
[rtg-users](https://groups.google.com/a/realtimegenomics.com/forum/#!forum/rtg-users)
discussion group is now available for general questions, tips, and
other discussions.

To be informed of new software releases, subscribe to the low-traffic
[rtg-announce](https://groups.google.com/a/realtimegenomics.com/forum/#!forum/rtg-announce)
group.

If you require a license for commercial use or wish to purchase
commercial support contact us via info@realtimegenomics.com.

---

## Prerequisites for building from source

* Java 1.8 or later
* apache ant 1.9 or later

## Check out source code for both RTG Tools and RTG Core

Building RTG Core requires the source code to both RTG Tools and RTG
Core, so you must clone both repositories:

    $ git clone https://github.com/RealTimeGenomics/rtg-tools.git
    $ git clone https://github.com/RealTimeGenomics/rtg-core.git
    $ cd rtg-core

To update, you will need to perform a `git pull` on both
repositories.  Advanced users may use `git subtree` to have both RTG
Tools and RTG Core within a single repository (it's what we use during
development).

## Compile / run unit tests

    $ ant runalltests

## Building RTG Core package

To build the RTG Core Non-Commercial package which can be locally
installed and run:

    $ ant zip-nojre

This will create an installation zip file under `dist`.

## Installation

Uncompress the installation zip:

    $ cd /my/install/dir/
    $ unzip /path/to/rtg-core/dist/rtg-core-VERSION-nojre.zip

Follow the instructions contained in the `README.txt`. This build will
use the system Java by default, so ensure it is Java 1.8 or later.

For a nice demonstration of the features of RTG Core for sex and
pedigree aware mapping and variant calling on data generated from
scratch using the RTG Core simulation tools, run the `demo-family.sh`
script contained in the scripts of the installation directory:

    $ cd /my/install/dir/rtg-core-VERSION/
    $ ./scripts/demo-family.sh $PWD/rtg

