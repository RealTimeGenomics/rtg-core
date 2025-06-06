/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.variant;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SingleMappedParams.SingleMappedParamsBuilder;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PosteriorUtils;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.DecomposerType;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;

import htsjdk.samtools.SAMFileHeader;

/**
 * A builder class for <code>VariantParams</code>.
 */
@TestClass("com.rtg.variant.VariantParamsTest")
public final class VariantParamsBuilder extends SingleMappedParamsBuilder<VariantParamsBuilder> {

  List<String> mImputedSamples = new ArrayList<>();
  Collection<File> mCalibrations = new ArrayList<>();
  int mQDefault = 20;
  int mMatedReadDefault = 20;
  int mUnmatedReadDefault = 20;
  int mMatedReadMax = 255;
  int mUnmatedReadMax = 255;
  boolean mIgnoreReadQualities = false;
  VariantOutputLevel mCallLevel = VariantOutputLevel.INTERESTING;
  double mMinAvrScore = 0;
  boolean mOutputNonSnps = true;
  GenomePriorParams mGenomePriors = null;
  double mInterestingThreshold = PosteriorUtils.unphredIfy(5.70342); // I.e. output uncertain ref calls where QUAL < ~5.7
  int mInterestingSeparation = 4;
  int mHyperComplexLength = 21;
  boolean mSimpleRepeatExtension = true;
  boolean mNoComplexCalls = false;
  CoverageThreshold mMaxCoverageFilter = new StaticThreshold(Integer.MAX_VALUE);
  CoverageThreshold mMaxCoverageBypass = new StaticThreshold(Integer.MAX_VALUE);
  CalibratedPerSequenceExpectedCoverage mExpectedCoverage = null;
  boolean mIgnoreQualityScores = false;
  //  boolean mCompleteGenomics = false;
  Double mMaxAmbiguity = null;
  Sex mSex = Sex.EITHER;
  ReferencePloidy mPloidy = ReferencePloidy.AUTO;
  String mMachineErrorName = null;
  boolean mVcfRp = false;
  boolean mOutputIndex = true;
  Calibrator mCalibrator = null;
  boolean mIonTorrent = false;
  int mChunkSize = 1000;
  int mLookAhead = 2;
  int mMaxReadLength = 1000;
  // ThreadingEnvironment mThreadingEnvironment = ThreadingEnvironment.SIMPLE;
  ThreadingEnvironment mThreadingEnvironment = ThreadingEnvironment.PARALLEL;
  Long mThreadingEnvironmentSeed = null;
  boolean mPruneHypotheses = false;
  double mIndelTriggerFraction = 0;
  DecomposerType mTrimSplit = DecomposerType.NONE;
  boolean mUsePropagatingPriors = false;

  boolean mExpandComplexReadQueries = false;

  boolean mComplexUseSoftClip = true;

  int mMaxComplexHypotheses = 6;

  int mMaxEmIterations = -1; // EmAlgorithm will turn this into DEFAULT_MAX_ITERATIONS
  File mPopulationPriorFile = null;

  File mAvrModelFile = null;

  GenomeRelationships mGenomeRelationships = null;
  GenomeConnectivity mGenomeConnectivity = null;
  double mNoDiseasePrior = 0.95;
  SAMFileHeader mUberHeader = null;
  ReferenceRanges<String> mReferenceRanges = null;
  File mRegionsFilterBedFile = null;

  double mMinVariantAllelicDepth = 0.0;
  double mMinVariantAllelicFraction = 0.0;

  EnumSet <VcfInfoField> mInfoAnnotations = EnumSet.noneOf(VcfInfoField.class);
  EnumSet<VcfFormatField> mFormatAnnotations = EnumSet.noneOf(VcfFormatField.class);
  SomaticParams mSomaticParams = new SomaticParamsBuilder().create();
  AlleleBalanceProbability mAlleleBalance = new NoAlleleBalance();
  RegionRestriction mForceComplexRegion;
  int mMinBaseQuality = 0;

  @Override
  protected VariantParamsBuilder self() {
    return this;
  }

  /**
   * Set the chunking size used internally for processing records.
   *
   * @param chunkSize chunk size
   * @return this builder
   */
  public VariantParamsBuilder chunkSize(final int chunkSize) {
    if (chunkSize < 1) {
      throw new IllegalArgumentException();
    }
    mChunkSize = chunkSize;
    return self();
  }

  /**
   * Set the number of chunks to pre-fetch.
   * @param lookahead number of chunks to fetch
   * @return this builder
   */
  public VariantParamsBuilder lookAhead(final int lookahead) {
    if (lookahead < 2) {
      throw new IllegalArgumentException();
    }
    mLookAhead = lookahead;
    return self();
  }

  /**
   * Set the maximum length for a read.
   *
   * @param maxReadLength the maximum length for a read.
   * @return this builder
   */
  public VariantParamsBuilder maxReadLength(final int maxReadLength) {
    if (maxReadLength < 1) {
      throw new IllegalArgumentException();
    }
    mMaxReadLength = maxReadLength;
    return self();
  }

  /**
   * Sets the names of samples to impute genotypes for.
   * @param samples the names of samples to impute.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder imputedSamples(final List<String> samples) {
    mImputedSamples = samples;
    return self();
  }

  /**
   * Sets the calibration input files.
   * @param calibrations the files containing calibrations.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder calibrations(final Collection<File> calibrations) {
    mCalibrations = calibrations;
    return self();
  }

  /**
   * Sets the default read mapping quality for mated reads, used when no
   * mapping quality is provided.
   * @param qual the phred scaled quality value. Default is 20.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder defaultMatedReadQuality(final int qual) {
    mMatedReadDefault = qual;
    return self();
  }

  /**
   * Sets the default read mapping quality for unmated reads, used when no
   * mapping quality is provided.
   * @param qual the phred scaled quality value. Default is 20.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder defaultUnmatedReadQuality(final int qual) {
    mUnmatedReadDefault = qual;
    return self();
  }

  /**
   * Sets the maximum read mapping quality for mated reads, used when no
   * mapping quality is provided.
   * @param qual the phred scaled quality value. Default is 255.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxMatedReadQuality(final int qual) {
    mMatedReadMax = qual;
    return self();
  }

  /**
   * Sets the default read mapping quality for unmated reads, used when no
   * mapping quality is provided.
   * @param qual the phred scaled quality value. Default is 20.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxUnmatedReadQuality(final int qual) {
    mUnmatedReadMax = qual;
    return self();
  }

  /**
   * @param ignore if true then ignore any supplied read quality scores (and just use the defaults).
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder ignoreReadQuality(final boolean ignore) {
    mIgnoreReadQualities = ignore;
    return self();
  }

  /**
   * Sets the default base quality, used when base quality values are not
   * provided.
   * @param qual the phred scaled quality value. Default is 20.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder defaultQuality(final int qual) {
    mQDefault = qual;
    return self();
  }

  /**
   * Sets whether a variant call is written at all positions.
   * @param flag true if calls should be written at all positions. Default is
   *          false.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder callLevel(final VariantOutputLevel flag) {
    mCallLevel = flag;
    return self();
  }
  /**
   * Sets whether a variant call is written at all positions.
   * @param score calls with AVR score below this will be failed
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder minAvrScore(final double score) {
    mMinAvrScore = score;
    return self();
  }

  /**
   * Sets whether indels should be called.
   * @param flag true if indels should be called. Default is false.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder outputNonSnps(final boolean flag) {
    mOutputNonSnps = flag;
    return self();
  }

  /**
   * Sets the priors probabilities to use.
   * @param priors builder object for priors.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder genomePriors(final GenomePriorParams priors) {
    mGenomePriors = priors;
    return self();
  }

  /**
   * A convenience method that takes the name of the priors resource file,
   * rather than a builder.
   * @param priors name of priors resource.
   * @return this builder, so calls can be chained.
   * @throws IOException if the priors file cannot be read.
   * @throws InvalidParamsException if priors file does not contain desired
   *           values.
   */
  public VariantParamsBuilder genomePriors(final String priors) throws IOException {
    mGenomePriors = GenomePriorParams.builder().genomePriors(priors).create();
    return self();
  }

  /**
   * Store the name of the over-ride machine errors.
   * @param val the name.
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder machineErrorName(final String val) {
    mMachineErrorName = val;
    return self();
  }

  /**
   * Sets the posterior below which an identity call is considered
   * <i>interesting</i>.
   * @param threshold posterior threshold for identity calls.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder interestingThreshold(final double threshold) {
    mInterestingThreshold = threshold;
    return self();
  }

  /**
   * Sets the maximum distance over which two <i>interesting</i> calls will be
   * considered part of the same complex region.
   * @param separation maximum distance
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder interestingSeparation(final int separation) {
    mInterestingSeparation = separation;
    return self();
  }

  /**
   * The length beyond which a complex region is considered Hyper Complex
   * @param length the length
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder hyperComplexLength(final int length) {
    mHyperComplexLength = length;
    return self();
  }

  /**
   * Whether to extend complex regions over simple repeats.
   * @param val whether to extend the regions.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder simpleRepeatExtension(final boolean val) {
    mSimpleRepeatExtension = val;
    return self();
  }

  /**
   * Whether to attempt calls in complex regions. If true then no such calls
   * are attempted but a single complex (type x) line is output.
   * @param val whether to attempt complex calls.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder noComplexCalls(final boolean val) {
    mNoComplexCalls = val;
    return self();
  }

  /**
   * @param threshold maximum coverage
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxCoverageFilter(final CoverageThreshold threshold) {
    mMaxCoverageFilter = threshold;
    return self();
  }

  /**
   * @param threshold maximum coverage for processing
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxCoverageBypass(final CoverageThreshold threshold) {
    mMaxCoverageBypass = threshold;
    return self();
  }

  /**
   * @param expectedCoverage the object containing the expected coverage information from the calibration files
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder expectedCoverage(final CalibratedPerSequenceExpectedCoverage expectedCoverage) {
    mExpectedCoverage = expectedCoverage;
    return self();
  }

  /**
   * Sets whether to ignore read quality scores and use the default instead.
   * @param flag true if read quality scores are to be ignored. Default is
   *          false.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder ignoreQualityScores(final boolean flag) {
    mIgnoreQualityScores = flag;
    return self();
  }

  /**
   * Sets the fraction of indels needed to trigger complex calling.
   * @param fraction the fraction of trivial indels needed to trigger complex calling.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder indelTriggerFraction(double fraction) {
    mIndelTriggerFraction = fraction;
    return self();
  }

  /**
   * Sets whether to use Ion Torrent specific processing.
   * @param flag true to use Ion Torrent specific processing.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder ionTorrent(final boolean flag) {
    mIonTorrent = flag;
    return self();
  }

  /**
   * @param flag true to use prune hypotheses during complex calling.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder pruneHypotheses(boolean flag) {
    mPruneHypotheses = flag;
    return self();
  }

  /**
   * Set threshold for ambiguity filtering.
   * @param th ratio to be used (if not set then apply no ambiguity filtering).
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxAmbiguity(final Double th) {
    mMaxAmbiguity = th;
    return self();
  }

  /**
   * Set the sex for the individual
   * @param sex male, female, or either (i.e. unknown)
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder sex(final Sex sex) {
    mSex = sex;
    return self();
  }

  /**
   * Set the default ploidy to use if no reference text file present.
   * @param ploidy diploid or haploid
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder ploidy(final ReferencePloidy ploidy) {
    mPloidy = ploidy;
    return self();
  }

  /**
   * Sets whether VCF output format should include the RTG posterior in the sample field
   * @param val the new value
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder vcfRp(final boolean val) {
    mVcfRp = val;
    return self();
  }

  /**
   * @param val true if TABIX index should be created for output SNP files.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder outputIndex(boolean val) {
    mOutputIndex = val;
    return self();
  }

  /**
   * @param c the calibrator if calibration has been supplied
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder calibrator(Calibrator c) {
    mCalibrator = c;
    return self();
  }

  /**
   * Create a new threading environment
   * @param threadingEnvironment threading environment
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder threadingEnvironment(ThreadingEnvironment threadingEnvironment) {
    mThreadingEnvironment = threadingEnvironment;
    return self();
  }

  /**
   * Create a new threading environment seed.
   * @param threadingEnvironmentSeed seed for random threading environment
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder threadingEnvironmentSeed(Long threadingEnvironmentSeed) {
    mThreadingEnvironmentSeed = threadingEnvironmentSeed;
    return self();
  }

  /**
   * @param value true if flag is enabled, false otherwise
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder trimSplit(DecomposerType value) {
    mTrimSplit = value;
    return self();
  }

  /**
   * The priors for no disease
   * @param p the priors
   * @return this, for chaining
   */
  public VariantParamsBuilder noDiseasePrior(final double p) {
    if (p <= 0 || p >= 1) {
      throw new IllegalArgumentException();
    }
    mNoDiseasePrior = p;
    return self();
  }

  /**
   * Set the SAM File uber header.
   * @param uberHeader the header that subsumes all input file headers
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder uberHeader(final SAMFileHeader uberHeader) {
    mUberHeader = uberHeader;
    return self();
  }

  /**
   * Set the reference ranges.
   * @param ranges the reference ranges that calling will operate over
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder referenceRanges(final ReferenceRanges<String> ranges) {
    mReferenceRanges = ranges;
    return self();
  }

  /**
   * Set the genome relationships.
   * @param genomeRelationships the genome relationships
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder genomeRelationships(final GenomeRelationships genomeRelationships) {
    mGenomeRelationships = genomeRelationships;
    return self();
  }

  /**
   * Set the population prior file
   * @param populationVcf the population prior VCF file
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder populationPriors(final File populationVcf) {
    mPopulationPriorFile = populationVcf;
    return self();
  }

  /**
   * Creates a VariantParams using the current builder configuration.
   * @return the new VariantParams
   */
  public VariantParams create() {
    return new VariantParams(self());
  }
  /**
   * @param val set to true to use forward backward algorithm
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder usePropagatingPriors(boolean val) {
    mUsePropagatingPriors = val;
    return self();
  }

  /**
   * Sets the maximum number of EM iterations.
   * @param iterations the number of EM iterations, or 0 to disable EM.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxEmIterations(int iterations) {
    mMaxEmIterations = iterations;
    return self();
  }

  /**
   * Sets the genome pedigree connectivity.
   * @param connectivity genome connectivity.
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder genomeConnectivity(GenomeConnectivity connectivity) {
    mGenomeConnectivity = connectivity;
    return self();
  }

  /**
   * Sets the hard cutoff on number of hypotheses
   * @param hypotheses the number of complex hypotheses to include
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder maxComplexHypotheses(int hypotheses) {
    mMaxComplexHypotheses = hypotheses;
    return self();
  }

  /**
   * @param value file containing an AVR model
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder avrModelFile(File value) {
    mAvrModelFile = value;
    return self();
  }

  /**
   * @param value the set of INFO annotations to output
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder infoAnnotations(EnumSet<VcfInfoField> value) {
    mInfoAnnotations = value;
    return self();
  }

  /**
   * @param value the set of FORMAT annotations to output
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder formatAnnotations(EnumSet<VcfFormatField> value) {
    mFormatAnnotations = value;
    return self();
  }

  /**
   * @param bedFile the bed file that will contain regions within which variants will have their FILTER field set
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder regionsFilterBedFile(File bedFile) {
    mRegionsFilterBedFile = bedFile;
    return self();
  }

  /**
   * @param depth output results with a variant allele depth meeting this threshold
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder minVariantAllelicDepth(double depth) {
    mMinVariantAllelicDepth = depth;
    return self();
  }

  /**
   * @param fraction output results with a variant allele fraction meeting this threshold
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder minVariantAllelicFraction(double fraction) {
    mMinVariantAllelicFraction = fraction;
    return self();
  }

  /**
   * @param params the configuration for a somatic caller
   * @return this builder, so calls can be chained.
   */
  public VariantParamsBuilder somaticParams(SomaticParams params) {
    mSomaticParams = params;
    return self();
  }

  /**
   * @param alleleBalance allele balance calculator
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder alleleBalance(AlleleBalanceProbability alleleBalance) {
    mAlleleBalance = alleleBalance;
    return self();
  }


  /**
   * @param val if true expand queries for reads by one base either side of a complex region
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder expandComplexReadQueries(boolean val) {
    mExpandComplexReadQueries = val;
    return self();
  }

  /**
   * @param val if true include soft clipped bases in complex calling evidence
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder complexUseSoftClip(boolean val) {
    mComplexUseSoftClip = val;
    return self();
  }

  /**
   * @param r region which will be forced to be a complex region during calling
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder forceComplexRegion(RegionRestriction r) {
    mForceComplexRegion = r;
    return self();
  }

  /**
   * @param regionString region which will be forced to be a complex region during calling. see {@link RegionRestriction#RegionRestriction(String)}
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder forceComplexRegion(String regionString) {
    mForceComplexRegion = new RegionRestriction(regionString);
    return self();
  }

  /**
   * @param minBaseQuality minimum base quality in phred scale
   * @return this builder, so calls can be chained
   */
  public VariantParamsBuilder minBaseQuality(int minBaseQuality) {
    mMinBaseQuality = minBaseQuality;
    return self();
  }
}
