#!/bin/bash -i

# This script runs a complete demonstration of RTG mapping and joint
# family variant calling, using data simulated on the fly.  When
# running it, give the full path to newly installed RTG, e.g.:
#
# demo-family.sh /path/to/rtg-core-NNN/rtg
#
# or (if rtg is installed on your $PATH)
#
# demo-family.sh rtg
#


if [ ! "$1" ]; then
    echo "Usage: $0 /path/to/rtg-core-NNN/rtg" >&2
    exit 1
fi

RTG=$1

# Define NOWAIT to run through without any pausing
# Define MD to format output as markdown

echo "Making directory for demo data: demo-family" >&2
if [ -e demo-family ]; then
    echo "Directory for working data 'demo-family' already exists, please delete it first" >&2
    exit 1
fi
if ! mkdir demo-family; then
    echo "Could not create a working directory for demo data.  Do you have write permission here?" >&2
    exit 1
fi
cd demo-family

echo "Checking RTG is executable" >&2
if ! "$RTG" version >/dev/null; then
    cat<<EOF >&2

Could not execute "$RTG" version. 

For this demo the path to RTG must be given as an absolute path.

EOF
    exit 1
fi


function pause() {
    if [ ! "$NOWAIT" ]; then
        echo
        read -ep "Press enter to continue..."
        echo
    fi
}

if [ "$MD" ]; then
    filter=("sed" "s/^/    /")
else
    filter=("cat")
fi

function docommand() {
    echo "\$ $@" | "${filter[@]}"
    "$@" 2>&1 | "${filter[@]}"
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        cat<<EOF >&2

Something unexpected happened during the demo. 
Check you are using the most recent version of RTG and this script. 
Exiting.
EOF
        exit 1
    fi
}

function doimage() {
    if [ "$MD" ]; then
        cat<<EOF
![Son1 ROC](file:///$PWD/$1)
EOF
    else
        cat<<EOF

You can look at the image 
$PWD/$1
with your favorite image viewer now or after the demo.

EOF
    fi
}

cat<<EOF

RTG Mapping and Variant Calling Demonstration
=============================================

In this demo we will give you a taste of the capabilities of RTG with
an end-to-end demonstration from simulated dataset generation through
mapping, variant calling, and analysis. 

To start with we will use RTG simulation utilities to generate a
synthetic dataset from scratch:

* genomesim - simulate a reference genome
* popsim - simulate population variants
* samplesim - generate two founder individuals
* childsim - simulate offspring of the two founders
* denovosim - simulate de novo mutations in some of the offspring
* readsim - simulate next-gen sequencing of the individuals

We will then use RTG alignment and variant calling commands to
demonstrate sex and pedigree-aware variant calling pipeline using the
following primary commands:

* map - read mapping to a reference
* snp - variant calling for a single sample
* family - joint variant calling for a single family

As well as some additional commands for quality assurance checks and
file manipulation:

* coverage - mapping coverage analysis
* chrstats - sample chromosome consistency checks
* mendelian - check variants for mendelian consistency
* vcffilter/vcfsubset - VCF record filtering
* vcfeval - compare two VCF call sets for agreement
* rocplot - produce static or interactive ROC graphs
* sdfstats - output information about data stored in SDF
* sdf2sam - extract reads from SDF in SAM

EOF
pause

cat<<EOF

Genome Simulation
-----------------

First we simulate a reference genome by generating random DNA, in this
case 10 chromosomes with lengths between 40kb and 50kb.  We will be
using fixed random number seeds during this demo in order to ensure we
have deterministic results.  (We take reproducability very seriously -
so you can be sure that you get repeatable results with RTG).

EOF
pause
docommand "$RTG" genomesim --output reference.sdf --num-contigs 10 --min-length 40000 --max-length 50000 --seed 42 --prefix Chr
cat<<EOF >reference.sdf/reference.txt
# Simulated reference. This file format is described in the user manual
version 1

# Default for any non-specified sequences, e.g. decoys
either  def     diploid linear

# Autosomes
either  seq     Chr1 diploid linear
either  seq     Chr2 diploid linear
either  seq     Chr3 diploid linear
either  seq     Chr4 diploid linear
either  seq     Chr5 diploid linear
either  seq     Chr6 diploid linear
either  seq     Chr7 diploid linear
either  seq     Chr8 diploid linear

# Sex chromosomes, here Chr9 is like human X, Chr10 is like human Y
male    seq     Chr9    haploid linear  Chr10
male    seq     Chr10   haploid linear  Chr9
female  seq     Chr9    diploid linear
female  seq     Chr10   none    linear

# If we were to define a pseudoautosomal region it might look like 
# this. While mapping and variant calling fully supports PAR regions,
# the simulation tools do not, so we won't use any in this demo.
#male    dup     Chr9:1501-4500  Chr10:5001-8000
EOF

pause
cat<<EOF

This command has created the reference as an SDF, which is a directory
structure that allows direct access to individual sequences or parts
of sequences, as well as containing metadata about the sequences
contained within.  RTG commands use SDFs for both reference and read
storage.

The reference SDF can optionally contain configuration that specifies
additional genomic information regarding ploidy and sex chromosomes,
we have also set that up, but here is what the configuration looks like:


EOF
docommand cat reference.sdf/reference.txt
pause

cat <<EOF

You can find out more information about any SDF by using the RTG
sdfstats command.  Here we will also request specific information
about how the reference sequences are interpreted for a male by adding
the flags '--sex male'.  Every RTG command accepts the --help option
in order to display the list of available options.

EOF
docommand "$RTG" sdfstats reference.sdf --sex male
pause



cat<<EOF

Recall that we have defined Chr9 as analogous to the human X
chromosome and Chr 10 as analogous to the human Y chromosome, and this
is reflected in the output.


Variant Simulation
------------------

Now we will simulate some variants on this reference, using the popsim
command.  The output will be a VCF containing these population
variants.

EOF
pause
docommand "$RTG" popsim --reference reference.sdf --seed 42 --output pop.vcf.gz
pause

cat<<EOF

The generated VCF contains variants with allele frequency
annotations that can be used to simulate a member of the population.
The types of variants and frequencies are based off variants from the
1000 Genomes Project.  We can examine some example variants using rtg
extract.  The extract command can be used with any SAM/BAM/BED/VCF
file that has been coordinate-sorted, block-compressed and indexed
according to standard NGS practise.  RTG commands automatically index
output files by default.

EOF
pause
docommand "$RTG" extract --header pop.vcf.gz Chr9:2000+1000
pause

cat<<EOF

Sample Simulation (including pedigree)
--------------------------------------

Now let's simulate a couple of members of the population that we can
use as parents for a family. For each sample we specify the desired
sample name and sex.  The samplesim command outputs a VCF including a
new sample column for the generated individual (and optionally an SDF
containing the whole genome for the sample, which we will use
later). We will run this twice, to generate each parent.

EOF
pause
docommand "$RTG" samplesim --reference reference.sdf --seed 572 --sex male --sample father --input pop.vcf.gz --output pop-1.vcf.gz --output-sdf genome-father.sdf
docommand "$RTG" samplesim --reference reference.sdf --seed 126 --sex female --sample mother --input pop-1.vcf.gz --output pop-2.vcf.gz --output-sdf genome-mother.sdf
pause

cat<<EOF

The genotypes selected for the two samples were determined on the
basis of the allele frequency annotations in the input VCF so a large
fraction of the low frequency population variants were not used.
Let's prune them from the VCF to keep things simple.

EOF
pause
docommand "$RTG" vcffilter --input pop-2.vcf.gz --output parents.vcf.gz --remove-all-same-as-ref
docommand "$RTG" extract --header parents.vcf.gz Chr9:2000+1000
pause

cat<<EOF

Now let's simulate some offspring of those two parents.  The RTG
childsim command obeys the reference chromosome information, selecting
variant genotypes following Mendelian inheritance and recombination.
For each child we specify the sample name of the father and mother,
and the desired sample name and sex of the child.

For the final child, we will use RTG denovosim to add novel variants
that are not present in either of the parents (the denovosim command
can also be used to simulate tumor/normal genomes).

As with the samplesim command, both childsim and denovosim generate
the individual as a new sample column in the output VCF, as well as
generating the whole genome SDF which we will use next for read
simulation.

EOF
pause
docommand "$RTG" childsim --reference reference.sdf --seed 837 --input parents.vcf.gz --output family-1.vcf.gz --output-sdf genome-son1.sdf --father father --mother mother --sex male --sample son1
docommand "$RTG" childsim --reference reference.sdf --seed 923 --input family-1.vcf.gz --output family-2.vcf.gz --output-sdf genome-son2.sdf --father father --mother mother --sex male --sample son2
docommand "$RTG" childsim --reference reference.sdf --seed 269 --input family-2.vcf.gz --output family-3.vcf.gz --output-sdf genome-daughter1.sdf --father father --mother mother --sex female --sample daughter1
docommand "$RTG" childsim --reference reference.sdf --seed 284 --input family-3.vcf.gz --output family-4.vcf.gz --father father --mother mother --sex female --sample daughter2-initial
docommand "$RTG" denovosim --reference reference.sdf --seed 841 --num-mutations 50 --input family-4.vcf.gz --output family-5.vcf.gz --output-sdf genome-daughter2.sdf --original daughter2-initial --sample daughter2
docommand "$RTG" vcfsubset --remove-sample daughter2-initial --input family-5.vcf.gz --output family.vcf.gz
rm pop-[0-9]*.vcf* family-[0-9]*.vcf*
pause

cat<<EOF

We can look at some of the variants produced for the samples on Chr9
(one of the sex chromosomes in our synthetic genome) and see that
those for the male samples are generated as haploid variants, as well
as the occasional de novo variant (annotated with a DN value of 'Y').

EOF
pause
docommand "$RTG" extract family.vcf.gz Chr9:10000+10000
pause

cat<<EOF

If we use rtg sdfstats --lengths to look at the genome SDF for one of
the male samples, we see that the chromosomes have been generated
according to the appropriate ploidy for the sample sex. Each
chromosome incorporates the appropriate variant alleles specified for
that sample in the VCF, so you will see slight variations in the
lengths of diploid pairs.

EOF
pause
docommand "$RTG" sdfstats --lengths genome-son1.sdf
pause



cat<<EOF

Read Simulation
---------------

Now we have a reference genome and have also generated simulated
genomes for the members of our family, we can simulate next generation
sequencing of the sample genomes, using RTG readsim.  We will simulate
2 x 100bp paired-end Illumina sequencing using a mean fragment size of
300bp, at a sequencing coverage of about 20x for most of our samples
(since our sample genome SDFs include two copies of each diploid
chromosome, we instruct the readsim command to use coverage 10).  One
of the samples will be generated at a lower coverage.  The simulated
reads include sequencer errors according to the selected machine
type, but in this case we will increase the error rates in order to
make the mapping and variant calling harder.  The results of read
simulation will be stored in an SDF containing the reads for each
sample.

EOF
pause
#readsim_opts="--machine=illumina_pe -L 100 -R 100 -m 200 -M 400 --coverage 10 --qual-range 2-20 --Xmnp-event-rate=0.02 --Xinsert-event-rate=0.005 --Xdelete-event-rate=0.005"
readsim_opts="--machine=illumina_pe -L 100 -R 100 -m 200 -M 400 --qual-range 2-20 --Xmnp-event-rate=0.02 --Xinsert-event-rate=0.0005 --Xdelete-event-rate=0.0005"
rgcommon="PL:ILLUMINA\tPI:300\tDS:Simulated dataset"
seed=5643
for genome in father mother son2 daughter1 daughter2; do
    seed=$[seed + 5]
    docommand "$RTG" readsim --input genome-$genome.sdf --output reads-$genome.sdf --seed $seed --sam-rg "@RG\tID:rg_$genome\tSM:$genome\t$rgcommon" --coverage 10 $readsim_opts|| exit 1
done
for genome in son1; do
    seed=$[seed + 5]
    docommand "$RTG" readsim --input genome-$genome.sdf --output reads-$genome.sdf --seed $seed --sam-rg "@RG\tID:rg_$genome\tSM:$genome\t$rgcommon" --coverage 5 $readsim_opts|| exit 1
done
pause

cat<<EOF

The RTG sdfstats command can also be used to retrieve information
about read sets.  We have included SAM read group information directly
in the SDF (this information could also be supplied at alignment
time).  Let's also extract some of the reads in SAM format using rtg
sdf2sam so you can see what they look like.

EOF

pause
docommand "$RTG" sdfstats reads-son1.sdf
docommand "$RTG" sdf2sam --input reads-son1.sdf --output - --start-id 0 --end-id 5
pause

cat<<EOF

Read Mapping
------------

Now we have set up our datasets we start the real work, of mapping the
reads to the reference genome.  When mapping a read set to the
reference, you should ensure that SAM read group information is
supplied, so that downstream variant calling can identify which sample
each alignment corresponds to.  When read data is being supplied from
an SDF or a SAM file that already contains read group information,
this is automatically used, otherwise it should be supplied as an
extra parameter on the command-line.

For the sex-aware pipeline you also need to indicate the sex of each
sample.  During mapping this can be specified by an explicit --sex
parameter, but it is easier when dealing with multiple samples to
store the sex and pedigree information in a PED file.  This is a
standard format, we will create one and show you now.

EOF
cat<<EOF >family.ped
1	father	0	0	1	0
1	mother	0	0	2	0
1	son1	father	mother	1	0
1	son2	father	mother	1	0
1	daughter1	father	mother	2	0
1	daughter2	father	mother	2	0
EOF
pause
docommand cat family.ped
pause

cat<<EOF

Typically when mapping real-world sequencing data, it will be provided
as FASTQ or BAM files, and the RTG map command can directly accept
both of these as input data.  However if the dataset is large enough
that it cannot be mapped within the available RAM of your compute
node, you will need to perform mapping in batches, using the
--start-read and --end-read parameters to the map command to select
which reads are being mapped.  If there will be more than a couple of
mapping batches required, it will be more efficient to convert the
read data to SDF using the RTG format command, as this avoids repeated
parsing of the input data, as well as being more memory-efficient.

During mapping, calibration information is computed regarding error
rates and typical genomic depth of coverage.  It is important that
when processing exome sequencing data, you supply the BED file
containing the target regions, so that this calibration information is
correctly computed.

There are many parameters you can supply to the map command to adjust
alignment accuracy and speed to suit your requirements.

Enough talk, let's map the read data.

EOF
pause
for genome in father mother son1 son2 daughter1 daughter2; do
    docommand "$RTG" map --template reference.sdf --input reads-$genome.sdf --output map-$genome --pedigree family.ped
done
pause

cat<<EOF

As you can see, at the end of each run a summary of mapping results is
displayed.  Inside the mapping output directory are the alignments, an
HTML report, and other informative files. Here is a quick example
graph of the fragment length distribution from one of the runs, which
allows you to check whether the minimum or maximum fragment length
parameters need adjusting (the full report contains this graph in
context along with other data).

EOF
doimage map-son1/index_files/Fragment-Length-Distribution.png

cat<<EOF

RTG also includes a coverage command which lets you analyze the depth
and breadth of coverage after mapping.  Let's run for each sample as a
sanity check. Like the map command, coverage command outputs a summary
in addition to the result files created in the output directory.

EOF
pause
for genome in father mother son1 son2 daughter1 daughter2; do
    docommand "$RTG" coverage --template reference.sdf --output coverage-$genome map-$genome/alignments.bam
done
pause

cat<<EOF

For example, the coverage command output includes the coverage
distribution and the cumulative coverage distribution so you can
examine the how much of each genome is covered by at least N depth of
coverage (these images are from the HTML report).

EOF
doimage coverage-son2/index_files/coverage.png
doimage coverage-son2/index_files/cumulative_coverage.png
pause

cat<<EOF

RTG also provides tools to assist with detecting sample mis-labelling
or incorrect pedigree information in two ways.

First, during mapping, chromosome coverage levels are checked for
consistency with each other given the expected levels for samples of
the specified sex.  A warning will be issued if inconsistencies are
detected.  This tool can also be invoked standalone, which is useful
when mapping has been split into several batches for each sample, or
to get a multi-sample summary.

EOF
pause
docommand "$RTG" chrstats --template reference.sdf --pedigree family.ped map-*/alignments.bam
pause

cat<<EOF

When a sample is detected as inconsistent, a warning is displayed
during mapping and the chrstats command may suggest what sex the
individual may actually be. Lets's show that now, by mapping one of
the samples incorrectly.

EOF
pause
docommand "$RTG" map --template reference.sdf --input reads-daughter1.sdf --output wrong-map-daughter1 --sex male
docommand "$RTG" chrstats --template reference.sdf --sex male wrong-map-daughter1/alignments.bam
pause

cat<<EOF

Examination of the output of the coverage command can be used to
obtain more detailed information about per-sequence coverage.  The
other consistency tool is used after variant calling to look for
abnormal mendelian inconsistencies or discordance between
samples. We'll come to that later.

EOF
pause

cat<<EOF

Variant Calling
---------------

RTG includes several commands for calling variants, depending on the
usage. The major callers are:

* snp - single-sample, human or non-human
* family - multiple samples within a single nuclear family. Multiple 
  children are permitted (it's not just for trios)
* population - multiple samples, which can contain related or unrelated 
  individuals according to pedigree.
* somatic - paired normal/tumor calling.

The multi-sample variant callers incorporate all samples jointly in
the Bayesian computations, which gives more accurate results than if
variant calling were carried out on each sample individually.

First we will run the single-sample caller on one of our samples as an
example and a basis for comparison. We just need to supply the
reference, pedigree file (to look up the sex of the sample), and the
BAM file created during mapping.

EOF
pause
docommand "$RTG" snp --template reference.sdf --pedigree family.ped --output snp-son1 map-son1/alignments.bam
pause

cat<<EOF

The variant caller provides summary statistics regarding the breakdown
of calls for each sample, as well as the primary VCF file and other
informative files within the output directory.

One thing that can be immediately noted from the summary statistics is
that this male sample includes haploid variant calls (made on the sex
chromosomes).

Now let's use the family caller.  In comparison with the single-sample
caller, we just need to supply the BAM files containing mappings of
the additional family members.

EOF
pause
docommand "$RTG" family --template reference.sdf --pedigree family.ped --output family-calls map-*/alignments.bam
pause

cat<<EOF

From these summary statistics you can see that child genotypes are
automatically phased in cases where this can be determined by
inheritance, and that putative de novo variant calls are also being
determined (and these are automatically marked as such in the output
VCF).

As a quick check that pedigree is being used to improve call
consistency, we will run the resulting VCF file through the RTG
mendelian command. This command outputs statistics on parental
concordance and calls inconsistent with Mendelian inheritance.

EOF
pause
docommand "$RTG" mendelian --template reference.sdf --input family-calls/snps.vcf.gz
pause

cat<<EOF

Variant Call Benchmarking and Comparison
----------------------------------------

Now let's compare the called variants with those in our gold standard,
the generated variants.  RTG includes a sophisticated variant
comparison tool called vcfeval, which is able to deal with the fact
that there are multiple equivalent ways to represent the same
variants, particularly those involving complex situations such as
block substitutions or indels.  Unlike other tools that compare
variant positions, alleles, and genotypes directly, RTG vcfeval
performs comparison at the level of haplotypes by replaying the
variants into the reference.  We will run vcfeval first to evaluate
the result of single-sample calling, and then to evaluate the family
called variants for the same sample.

EOF
pause
docommand "$RTG" vcfeval --template reference.sdf --baseline family.vcf.gz --calls snp-son1/snps.vcf.gz --sample son1 -o vcfeval-son1-snp
docommand "$RTG" vcfeval --template reference.sdf --baseline family.vcf.gz --calls family-calls/snps.vcf.gz --sample son1 -o vcfeval-son1-family
pause

cat<<EOF

We can see from the summary reports that the family caller was able to
make use of the other family members to increase the calling accuracy.

As well as a summary of the overall calling accuracy, vcfeval creates
an output file for performing ROC analysis.  Variant callers typically
err on the side of outputting too many variants (favoring
sensitivity), since it's always possible to post-filter to remove
false positives.  ROC analysis assists with choosing appropriate
filter thresholds, by demonstrating the trade-off between sensitivity
and precision with respect to a filter threshold.  The default is to
analyse using the FORMAT GQ field, but when you run vcfeval, you can
specify the --vcf-score-field option to use a different attribute
(such as QUAL or AVR)

The RTG rocplot command allows ROC plots for multiple ROC data files
to be viewed using either an interactive gui, or via creation of a
static image.  We'll generate a static image now, (As an alternative,
you could also decompress weighted_roc.tsv.gz file and load it as
tab-separated file into a spreadsheet or other graphing package.)

EOF
pause
docommand "$RTG" rocplot --title ROC-Son1 vcfeval-son1-*/weighted_roc.tsv.gz --png=roc-son1.png
doimage roc-son1.png

cat<<EOF

That one isn't too interesting, since this is simulated data with a
relatively low number of data points, so the ROC tracks the top-left
pretty well.  The de novo calling is much trickier, since a de novo
variant only has direct evidence from the single sample containing
that variant, so let's see how we did there.

EOF
pause
docommand "$RTG" vcffilter --input family.vcf.gz --output family-denovo.vcf.gz --min-denovo-score=0 --sample daughter2 
docommand "$RTG" vcfeval --template reference.sdf --baseline family-denovo.vcf.gz --calls family-calls/snps.vcf.gz --sample daughter2 --output vcfeval-daughter2-denovo --vcf-score-field DNP
pause

cat<<EOF

Of all the variants that were suggested to be de novo, there were 0
false positive and only a couple of false negatives.  We included the
full call-set in the comparison, so there were a lot of calls that
would be false positives with respect to the list of de novos that we
generated, so focusing effort on those with good DNP scores makes
sense when hunting for those de novo variants.

That's it for the demo.  Feel free to look around in the various
output files that were created inside this demo-family directory, and
try out RTG on some real data of your own.

EOF


