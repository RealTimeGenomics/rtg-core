/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SingleMappedParams;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.test.params.ParamsNoField;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;

import htsjdk.samtools.SAMFileHeader;

/**
 */
public final class VariantParams extends SingleMappedParams implements VariantOutputOptions, Integrity {

  /** Name of file used for output of called SNPs in VCF format. */
  public static final String VCF_OUT_SUFFIX = "snps.vcf";
  /** Name of file used for output of complex regions in BED format. */
  public static final String BED_OUT_SUFFIX = "regions.bed";

  /**
   * Creates a VariantParams builder.
   * @return the builder.
   */
  public static VariantParamsBuilder builder() {
    return new VariantParamsBuilder();
  }

  private final int mQDefault;
  private final int mMatedReadDefault;
  private final int mUnmatedReadDefault;
  private final int mMatedReadMax;
  private final int mUnmatedReadMax;
  private final boolean mIgnoreReadQualities;
  private final List<String> mImputedSamples;
  private final Collection<File> mCalibrations;
  private final VariantOutputLevel mCallLevel;
  private final boolean mOutputNonSnps;
  private final GenomePriorParams mGenomePriors;
  private final double mInterestingThreshold;
  private final int mIndelTriggerFraction;
  private final int mInterestingSeparation;
  private final int mHyperComplexLength;
  private final boolean mSimpleRepeatExtension;
  private final boolean mNoComplexCalls;
  private final CoverageThreshold mMaxCoverageFilter;
  private final CoverageThreshold mMaxCoverageBypass;
  private final CalibratedPerSequenceExpectedCoverage mExpectedCoverage;
  private final boolean mIgnoreQualityScores;
  private final Double mMaxAmbiguity;
  private final Sex mSex;
  private final ReferencePloidy mPloidy;
  private final String mMachineErrorName;
  private final boolean mVcfRp;
  private final boolean mOutputIndex;
  private final Calibrator mCalibrator;
  private final int mChunkSize;
  private final int mLookAhead;
  private final int mMaxReadLength;
  private final boolean mIonTorrent;
  private final ThreadingEnvironment mThreadingEnvironment;
  private final Long mThreadingEnvironmentSeed;
  private final boolean mPruneHypotheses;
  private final boolean mEnableTrimSplit;
  private final File mPopulationPriorFile;
  private final int mMaxEmIterations;
  private final int mMaxComplexHypotheses;
  private final File mAvrModelFile;
  private final double mMinAvrScore;
  private final SAMFileHeader mUberHeader;
  private final ReferenceRanges<String> mReferenceRanges;
  private final GenomeRelationships mGenomeRelationships;
  private final GenomeConnectivity mGenomeConnectivity;
  private final double mNoDiseasePrior;
  private final boolean mUsePropagatingPriors;
  private final EnumSet<VcfInfoField> mInfoAnnotations;
  private final EnumSet<VcfFormatField> mFormatAnnotations;
  private final File mRegionsFilterBedFile;
  private final int mMinVariantAlleleCount;
  private final double mMinVariantAlleleFraction;
  private final SomaticParams mSomaticParams;
  private final AlleleBalanceProbability mAlleleBalance;
  private final boolean mExpandComplexReadQueries;

  /**
   * @param builder the builder object.
   */
  VariantParams(final VariantParamsBuilder builder) {
    super(builder);
    mQDefault = builder.mQDefault;
    mMatedReadDefault = builder.mMatedReadDefault;
    mUnmatedReadDefault = builder.mUnmatedReadDefault;
    mMatedReadMax = builder.mMatedReadMax;
    mUnmatedReadMax = builder.mUnmatedReadMax;
    mIgnoreReadQualities = builder.mIgnoreReadQualities;
    mImputedSamples = builder.mImputedSamples;
    mCalibrations = builder.mCalibrations;
    mCallLevel = builder.mCallLevel;
    mOutputNonSnps = builder.mOutputNonSnps;
    mGenomePriors = builder.mGenomePriors;
    mInterestingThreshold = builder.mInterestingThreshold;
    mInterestingSeparation = builder.mInterestingSeparation;
    mHyperComplexLength = builder.mHyperComplexLength;
    mSimpleRepeatExtension = builder.mSimpleRepeatExtension;
    mNoComplexCalls = builder.mNoComplexCalls;
    mMaxCoverageFilter = builder.mMaxCoverageFilter;
    mMaxCoverageBypass = builder.mMaxCoverageBypass;
    mExpectedCoverage = builder.mExpectedCoverage;
    mIgnoreQualityScores = builder.mIgnoreQualityScores;
    mMachineErrorName = builder.mMachineErrorName;
    mMaxAmbiguity = builder.mMaxAmbiguity;
    mSex = builder.mSex;
    mPloidy = builder.mPloidy;
    mVcfRp = builder.mVcfRp;
    mOutputIndex = builder.mOutputIndex;
    mCalibrator = builder.mCalibrator;
    mChunkSize = builder.mChunkSize;
    mLookAhead = builder.mLookAhead;
    mMaxReadLength = builder.mMaxReadLength;
    mIonTorrent = builder.mIonTorrent;
    mThreadingEnvironment = builder.mThreadingEnvironment;
    mThreadingEnvironmentSeed = builder.mThreadingEnvironmentSeed;
    mPruneHypotheses = builder.mPruneHypotheses;
    mIndelTriggerFraction = builder.mIndelTriggerFraction;
    mEnableTrimSplit = builder.mEnableTrimSplit;
    mAvrModelFile = builder.mAvrModelFile;
    mMinAvrScore = builder.mMinAvrScore;
    mUberHeader = builder.mUberHeader;
    mReferenceRanges = builder.mReferenceRanges;
    mGenomeRelationships = builder.mGenomeRelationships;
    mGenomeConnectivity = builder.mGenomeConnectivity;
    mNoDiseasePrior = builder.mNoDiseasePrior;
    mPopulationPriorFile = builder.mPopulationPriorFile;
    mUsePropagatingPriors = builder.mUsePropagatingPriors;
    mMaxEmIterations = builder.mMaxEmIterations;
    mMaxComplexHypotheses = builder.mMaxComplexHypotheses;
    mInfoAnnotations = builder.mInfoAnnotations;
    mFormatAnnotations = builder.mFormatAnnotations;
    mRegionsFilterBedFile = builder.mRegionsFilterBedFile;
    mMinVariantAlleleCount = builder.mMinVariantAlleleCount;
    mMinVariantAlleleFraction = builder.mMinVariantAlleleFraction;
    mSomaticParams = builder.mSomaticParams;
    mAlleleBalance = builder.mAlleleBalance;
    mExpandComplexReadQueries = builder.mExpandComplexReadQueries;
  }

  @Override
  public VariantOutputLevel callLevel() {
    return mCallLevel;
  }

  /**
   * Check if non-SNPs should be included in output.
   * @return true iff non-SNPs should be included in output.
   */
  public boolean outputNonSnps() {
    return mOutputNonSnps;
  }

  /**
   * Get default quality value for individual nucleotides.
   * @return default quality value.
   */
  public int qDefault() {
    return mQDefault;
  }

  /**
   * @return Default quality for mated reads.
   */
  public int matedReadDefault() {
    return mMatedReadDefault;
  }

  /**
   * @return Default quality for unmated reads.
   */
  public int unmatedReadDefault() {
    return mUnmatedReadDefault;
  }

  /**
   * @return Default quality for mated reads.
   */
  public int matedReadMax() {
    return mMatedReadMax;
  }

  /**
   * @return Default quality for unmated reads.
   */
  public int unmatedReadMax() {
    return mUnmatedReadMax;
  }

  /**
   * @return true if read qualities to be ignored.
   */
  public boolean ignoreReadQualities() {
    return mIgnoreReadQualities;
  }

  /**
   * @return file used for stream returned by {@link #vcfStream()}
   */
  @ParamsNoField
  public File vcfFile() {
    return outFile(VCF_OUT_SUFFIX);
  }

  /**
   * @return file used for stream returned by {@link #bedStream()}
   */
  @ParamsNoField
  public File bedFile() {
    return outFile(BED_OUT_SUFFIX);
  }

  /**
   * @return the stream for writing the VCF file.
   * @throws IOException whenever.
   */
  @ParamsNoField
  public OutputStream vcfStream() throws IOException {
    return outStream(VCF_OUT_SUFFIX);
  }

  /**
   * @return the stream for writing the bed file of complex regions.
   * @throws IOException whenever.
   */
  @ParamsNoField
  public OutputStream bedStream() throws IOException {
    return outStream(BED_OUT_SUFFIX);
  }

  /**
   * Get information about all prior probabilities.
   * @return A GenomePriorParams object.
   */
  public GenomePriorParams genomePriors() {
    return mGenomePriors;
  }

  /**
   * @return the name of the over-ride machine errors.
   */
  public String machineErrorName() {
    return mMachineErrorName;
  }

  /**
   * Return the posterior threshold below which an identity call is considered
   * <i>interesting</i>.
   * @return posterior threshold.
   */
  public double interestingThreshold() {
    return mInterestingThreshold;
  }

  /**
   * Return the maximum distance over which two <i>interesting</i> calls will be
   * considered part of the same complex region.
   * @return interesting separation distance.
   */
  public int interestingSeparation() {
    return mInterestingSeparation;
  }

  /**
   * @return the length beyond which a complex region is considered Hyper Complex
   */
  public int hyperComplexLength() {
    return mHyperComplexLength;
  }

  /** @return true if complex regions should be extended over simple repeats */
  public boolean simpleRepeatExtension() {
    return mSimpleRepeatExtension;
  }

  @Override
  public boolean noComplexCalls() {
    return mNoComplexCalls;
  }

  /**
   * Return the maximum coverage for filtering output.
   * @return maximum read coverage
   */
  @Override
  public CoverageThreshold maxCoverageFilter() {
    return mMaxCoverageFilter;
  }

  /**
   * Return the maximum coverage to attempt processing.
   * @return maximum read coverage
   */
  public CoverageThreshold maxCoverageBypass() {
    return mMaxCoverageBypass;
  }

  /**
   * Return the expected coverage information.
   * @return expected read coverage information.
   */
  public CalibratedPerSequenceExpectedCoverage expectedCoverage() {
    return mExpectedCoverage;
  }

  /**
   * Return whether to ignore read quality scores, and use default quality score
   * instead.
   * @return ignore read quality scores
   */
  public boolean ignoreQualityScores() {
    return mIgnoreQualityScores;
  }

  /**
   * Return whether to output the posterior for the total variants as well as
   * the best variant.
   * @return true iff the total variant posterior is to be output.
   */
  @ParamsNoField
  public boolean nonidentityPosterior() {
    return true; //always get the non-identity posterior if we're outputting VCF
  }

  @Override
  public Double maxAmbiguity() {
    return mMaxAmbiguity;
  }

  /**
   * Return the sex of this individual.
   * @return the sex
   */
  public Sex sex() {
    return mSex;
  }

  /**
   * Return the default ploidy to use if no reference text file present.
   * @return the default ploidy
   */
  public ReferencePloidy ploidy() {
    return mPloidy;
  }

  /**
   * @return true if VCF output should include the RTG posterior in the sample field
   */
  public boolean vcfRp() {
    return mVcfRp;
  }

  /**
   * Return the names of samples to be imputed.
   * @return the names of imputed samples
   */
  public List<String> imputedSamples() {
    return mImputedSamples;
  }

  /**
   * Return the calibrations.
   * @return the calibrations
   */
  public Collection<File> calibrations() {
    return mCalibrations;
  }

  /**
   * @return true if TABIX index should be output.
   */
  public boolean outputIndex() {
    return mOutputIndex;
  }

  /**
   * @return the calibrator generated from given calibration files
   */
  public Calibrator calibrator() {
    return mCalibrator;
  }

  /** @return the fraction of indels needed to trigger complex calling */
  public int indelTriggerFraction() {
    return mIndelTriggerFraction;
  }

  /** @return true if ion torrent data is detected during this run */
  public boolean ionTorrent() {
    return mIonTorrent;
  }

  /** @return true if pruning hypotheses is desired */
  public boolean pruneHypotheses() {
    return mPruneHypotheses;
  }

  /** @return threading environment for this run */
  public ThreadingEnvironment threadingEnvironment() {
    return mThreadingEnvironment;
  }

  /** @return threading environment aware look ahead */
  public int threadingLookAhead() {
    final int lookAhead;
    switch (threadingEnvironment()) {
      case SINGLE:
        lookAhead = 1;
        break;
      case RANDOM:
      case PARALLEL:
        lookAhead = lookAhead() * execThreads() + 1; //The LOOKAHEAD_MULTIPLIER and 1 are heuristic
        break;
      default:
        throw new RuntimeException();
    }
    return lookAhead;
  }

  /** @return seed for threading environment, can be null */
  public Long threadingEnvironmentSeed() {
    return mThreadingEnvironmentSeed;
  }

  /**
   * @return the number of chunks to prefetch.
   */
  public int lookAhead() {
    return mLookAhead;
  }

  /**
   * @return the number of positions to include when dividing work into chunks.
   */
  public int chunkSize() {
    return mChunkSize;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mapped());
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    Exam.assertTrue(integrity() && (mGenomePriors == null || mGenomePriors.globalIntegrity()));
    return true;
  }

  /**
   * @return the maximum length of a read.
   */
  public int maxReadLength() {
    return mMaxReadLength;
  }

  /**
   * @return enable trimming and splitting
   */
  public boolean enableTrimSplit() {
    return mEnableTrimSplit;
  }

  /**
   * @return get the SAM uber header (an aggregate form of the individual sam file headers)
   */
  public SAMFileHeader uberHeader() {
    return mUberHeader;
  }

  /**
   * @return the reference ranges that calling will operate over
   */
  public ReferenceRanges<String> referenceRanges() {
    return mReferenceRanges;
  }

  /**
   * @return get the genome relationships
   */
  public GenomeRelationships genomeRelationships() {
    return mGenomeRelationships;
  }

  /**
   * @return return the genome pedigree connectivity
   */
  public GenomeConnectivity genomeConnectivity() {
    return mGenomeConnectivity;
  }

  /**
   * @return prior probability that a position does not explain a disease.
   */
  public double noDiseasePrior() {
    return mNoDiseasePrior;
  }

  /**
   * @return VCF file containing population priors
   */
  public File populationPriorFile() {
    return mPopulationPriorFile;
  }

  /**
   * @return the maximum number of EM iterations to perform. 0 indicates EM is disabled
   */
  public int maxEmIterations() {
    return mMaxEmIterations;
  }

  /**
   * @return the maximum number of complex hypotheses to consider
   */
  public int maxComplexHypotheses() {
    return mMaxComplexHypotheses;
  }

  /**
   * @return true if should use forward backward algorithm
   */
  public boolean usePropagatingPriors() {
    return mUsePropagatingPriors;
  }

  /**
   * @return file for AVR model to apply, null if no file set
   */
  public File avrModelFile() {
    return mAvrModelFile;
  }

  /**
   * Fail calls below this score
   * @return a double in the AVR score range [0-1]
   */
  public double minAvrScore() {
    return mMinAvrScore;
  }

  /**
   * @return set of INFO annotations to output
   */
  public EnumSet<VcfInfoField> infoAnnotations() {
    return mInfoAnnotations;
  }

  /**
   * @return set of FORMAT annotations to output
   */
  public EnumSet<VcfFormatField> formatAnnotations() {
    return mFormatAnnotations;
  }

  /**
   * @return bed file containing restriction regions used to set the FILTER status of variants. (not used to influence calling regions).
   */
  public File regionsFilterBedFile() {
    return mRegionsFilterBedFile;
  }

  /**
   * @return the minimum variant allelic count to output a call. Calls below this threshold will be filtered
   */
  public int minVariantAlleleCount() {
    return mMinVariantAlleleCount;
  }

  /**
   * @return the minimum variant allelic fraction output a call. Calls below this threshold will be filtered
   */
  public double minVariantAlleleFraction() {
    return mMinVariantAlleleFraction;
  }
  /**
   * @return the params for somatic specific callers
   */
  public SomaticParams somaticParams() {
    return mSomaticParams;
  }
  /**
   * @return the allele balance probability calculator
   */
  public AlleleBalanceProbability alleleBalance() {
    return mAlleleBalance;
  }

  /**
   * @return if true expand queries for reads by one base either side of a complex region
   */
  public boolean expandComplexReadQueries() {
    return mExpandComplexReadQueries;
  }

  /**
   * Create a builder with all the values set to those of this object.
   * @return a builder
   */
  public VariantParamsBuilder cloneBuilder() {
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    cloneSet(vpb);
    return vpb;
  }

  private void cloneSet(final VariantParamsBuilder vpb) {
    vpb
    .name(name())
    .chunkSize(mChunkSize)
    .maxReadLength(mMaxReadLength)
    .calibrator(mCalibrator)
    .outputIndex(mOutputIndex)
    .vcfRp(mVcfRp)
    .machineErrorName(mMachineErrorName)
    .sex(sex())
    .maxAmbiguity(mMaxAmbiguity)
    .ignoreQualityScores(mIgnoreQualityScores)
    .maxCoverageFilter(mMaxCoverageFilter)
    .maxCoverageBypass(mMaxCoverageBypass)
    .expectedCoverage(mExpectedCoverage)
    .noComplexCalls(mNoComplexCalls)
    .hyperComplexLength(mHyperComplexLength)
    .interestingSeparation(mInterestingSeparation)
    .interestingThreshold(mInterestingThreshold)
    .genomePriors(mGenomePriors)
    .outputNonSnps(mOutputNonSnps)
    .callLevel(mCallLevel)
    .filterParams(filterParams())
    .imputedSamples(mImputedSamples)
    .calibrations(mCalibrations)
    .ignoreReadQuality(mIgnoreReadQualities)
    .outputParams(outputParams())
    .defaultQuality(mQDefault)
    .defaultMatedReadQuality(mMatedReadDefault)
    .defaultUnmatedReadQuality(mUnmatedReadDefault)
    .maxMatedReadQuality(mMatedReadMax)
    .maxUnmatedReadQuality(mUnmatedReadMax)
    .pruneHypotheses(mPruneHypotheses)
    .threadingEnvironment(mThreadingEnvironment)
    .threadingEnvironmentSeed(mThreadingEnvironmentSeed)
    .ploidy(mPloidy)
    .enableTrimSplit(mEnableTrimSplit)
    .genomeRelationships(mGenomeRelationships)
    .noDiseasePrior(mNoDiseasePrior)
    .ionTorrent(mIonTorrent)
    .indelTriggerFraction(mIndelTriggerFraction)
    .execThreads(execThreads())
    .ioThreads(ioThreads())
    .mapped(mapped())
    .uberHeader(uberHeader())
    .genome(genome())
    .populationPriors(mPopulationPriorFile)
    .usePropagatingPriors(mUsePropagatingPriors)
    .maxEmIterations(mMaxEmIterations)
    .maxComplexHypotheses(mMaxComplexHypotheses)
    .avrModelFile(mAvrModelFile)
    .minAvrScore(mMinAvrScore)
    .infoAnnotations(EnumSet.copyOf(mInfoAnnotations))
    .formatAnnotations(EnumSet.copyOf(mFormatAnnotations))
    .regionsFilterBedFile(mRegionsFilterBedFile)
    .referenceRanges(mReferenceRanges)
    .minVariantAlleleCount(mMinVariantAlleleCount)
    .minVariantAlleleFraction(mMinVariantAlleleFraction)
    .somaticParams(somaticParams())
    .alleleBalance(alleleBalance())
    ;
  }

  @Override
  public String toString() {
    final String pref = "    ";
    final StringBuilder sb = new StringBuilder();
    sb.append("VariantParams");
    final Collection<File> mapped = mapped();
    if (mapped != null) {
      sb.append(" mapped reads=[");
      int i = 0;
      for (final File file : mapped) {
        if (i++ > 0) {
          sb.append(", ");
        }
        sb.append(file.getPath());
      }
      sb.append("]");
    }
    sb.append(" q_default=").append(mQDefault)
      .append(" mated_read_default=").append(mMatedReadDefault)
      .append(" unmated_read_default=").append(mUnmatedReadDefault).append(LS);
    sb.append(" mated_read_max=").append(mMatedReadMax)
      .append(" unmated_read_max=").append(mUnmatedReadMax)
      .append(" ignore_map_qualities=").append(mIgnoreReadQualities).append(LS);
    sb.append(" hypercomplex_length=").append(mHyperComplexLength)
      .append(" non_identity_posterior=").append(nonidentityPosterior()).append(LS);
    sb.append(" machine=").append(mMachineErrorName)
      .append(" vcf_rp=").append(mVcfRp)
      .append(" output_index=").append(mOutputIndex).append(LS);
    sb.append(" call_level=").append(mCallLevel)
      .append(" indels=").append(mOutputNonSnps)
      .append(" noComplexCalls=").append(mNoComplexCalls).append(LS);
    sb.append(" interesting_threshold=").append(Utils.realFormat(mInterestingThreshold, 4)).append(LS);
    sb.append(" interesting_separation=").append(mInterestingSeparation).append(LS);
    sb.append(" indel_trigger_fraction=").append(mIndelTriggerFraction).append(LS)
      .append(" max_coverage_filter=").append(mMaxCoverageFilter).append(LS)
      .append(" max_coverage_bypass=").append(mMaxCoverageBypass).append(LS);
    sb.append(" ignore_quality_scores=").append(mIgnoreQualityScores).append(LS);
    sb.append(" max_ambiguity=").append(mMaxAmbiguity).append(LS);
    sb.append(" sex=").append(mSex).append(LS);
    sb.append(" ploidy=").append(mPloidy).append(LS);
    sb.append(" chunk_size=").append(mChunkSize)
      .append(" lookahead=").append(mLookAhead)
      .append(" max_read_length=").append(mMaxReadLength).append(LS);
    sb.append(" threading_environment=").append(mThreadingEnvironment)
      .append(" treading_environment_seed=").append(mThreadingEnvironmentSeed).append(LS);
    sb.append(" exec_threads=").append(execThreads())
      .append(" io_threads=").append(ioThreads()).append(LS);
    sb.append(" hyper_complex_threshold=").append(mHyperComplexLength).append(LS);
    sb.append(" ionTorrent=").append(mIonTorrent)
      .append(" prune_hypothesis=").append(mPruneHypotheses)
      .append(" enable_trim_split=").append(mEnableTrimSplit).append(LS);  //note currently writing IonTorrent out here is pointless as it is modified after this is printed.
    sb.append(" noDiseasePrior=").append(noDiseasePrior()).append(StringUtils.LS);
    sb.append(" Relationships:").append(genomeRelationships()).append(LS);
    sb.append(" max_em_iterations=").append(mMaxEmIterations).append(LS);
    sb.append(" genome_connectivity=").append(mGenomeConnectivity == null ? "null" : mGenomeConnectivity.toString()).append(LS);
    if (mPopulationPriorFile != null) {
      sb.append(" Population prior=").append(mPopulationPriorFile.getPath()).append(LS);
    }
    sb.append(mGenomePriors == null ? "null" : mGenomePriors.toString()).append(LS);
    if (genome() != null) {
      sb.append(pref).append(genome().toString()).append(LS);
    }
    if (outputParams() != null) {
      sb.append(pref).append(outputParams().toString()).append(LS);
    }
    if (filterParams() != null) {
      sb.append(pref).append(filterParams().toString()).append(LS);
    }
    sb.append(" AVR model: ").append(mAvrModelFile == null ? "Not set" : mAvrModelFile.getPath()).append(LS);
    sb.append(" min AVR score=").append(mMinAvrScore).append(LS);
    sb.append(" max_complex_hypotheses=").append(mMaxComplexHypotheses).append(LS);
    sb.append(" regions_bed_file=").append(mRegionsFilterBedFile).append(LS);
    sb.append(" min_variant_allele_count=").append(mMinVariantAlleleCount).append(LS);
    sb.append(" min_variant_allele_fraction=").append(mMinVariantAlleleFraction).append(LS);
    sb.append(" somatic_params=").append(mSomaticParams).append(LS);
    sb.append(" allele_balance=").append(mAlleleBalance).append(LS);
    return sb.toString();
  }
}
