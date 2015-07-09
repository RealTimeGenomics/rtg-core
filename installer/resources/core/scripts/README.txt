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


RTG Example Reference Files
===========================

These files illustrate how to specify a reference.txt file that is
used inside a formatted reference to allow sex-aware mapping, variant
calling, and analysis. For standard human references choose one of
these examples that matches the version and chromosome naming
convention of your reference.  You may need to manually adjust the
chromosome names if you are using NCBI identifiers rather than chrNN
or NN style. See the user manual for more information on reference.txt
specification and how to carry out sex-aware mapping and variant
calling.

hg18.example.reference.txt
hg19.example.reference.txt
GRCh37.example.reference.txt
GRCh38.example.reference.txt


Additional Scripts
==================

This directory also contains additional helper and demonstration
scripts. Each script contains a header with brief documentation as to
the purpose and usage.



