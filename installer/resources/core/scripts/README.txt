RTG Example Scripts
===================

Several subdirectories in this folder contain example scripts for
metagenomic analysis that are intended to be modified to suit your
specific pipeline needs.

sample-clustering        Visualization of sample clustering of 
                         similarity output using R.
species-abundance        Visualization of species output using Excel.
species-comparison       Visualization of multiple species outputs 
                         using R.

A README.txt in each directory contains more information.


RTG Example Reference File
==========================

These files illustrate how to specify a reference.txt file to allow
sex-specific mapping and variant calling. You may need to adjust the
chromosome names if you are using NCBI identifiers rather than chrNN
style. See the user manual for more information on reference.txt
specification and how to carry out sex-specific mapping and variant
calling.

hg18.example.reference.txt        An example reference.txt for hg18.
hg19.example.reference.txt        An example reference.txt for hg19.
