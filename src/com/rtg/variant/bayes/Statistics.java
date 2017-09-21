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

package com.rtg.variant.bayes;

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;

/**
 * Maintains various counts etc. that are used by the model logic during output.
 * It is expected that there will be classes of this that maintain other information that
 * is presented in human readable form.
 */
@TestClass("com.rtg.variant.bayes.snp.StatisticsSnpTest")
public abstract class Statistics<T extends AlleleStatistics<T>> implements Cloneable {

  protected T mCounts;
  protected Description mDescription;

  // This is used for short circuit triggering. Integer is OK
  private double mNonRefCount;

  // These stats are currently using integer counts, but are at least self-consistent
  // They could probably be pushed into specialized double / int versions if we wanted.
  private int mCountLeft;
  private int mCountRight;
  private int mCountUnmapped;

  // Workspace for storing allele indexes (avoids churn)
  private final int[] mAlleleIndexes = new int[2];

  /**
   * @param description about which statistics are being collected.
   * @param alleleStats holder for per-allele statistics
   */
  public Statistics(Description description, T alleleStats) {
    mDescription = description;
    mCounts = alleleStats;
  }

  /**
   * Increment the statistics.
   * @param evidence with probabilities for hypotheses.
   * @param reference code for reference allele
   */
  public final void increment(final EvidenceInterface evidence, int reference) {
    if (evidence.isUnmapped()) {
      ++mCountUnmapped;
      return; // Don't increment any other stats for unmapped evidence
    }

    final int read = evidence.read();
    if (read == Hypotheses.NO_HYPOTHESIS) {
      return;
    }

    // Used for short circuit detection, integer increments OK
    if (read != reference) {
      ++mNonRefCount;
    }

    // Used for internally self-consistent stats calculations, so can be integer increments
    final int readBasesLeft = evidence.getReadBasesLeft();
    final int readBasesRight = evidence.getReadBasesRight();
    final int left = readBasesLeft >= readBasesRight ? 1 : 0;
    final int right = readBasesRight >= readBasesLeft ? 1 : 0;
    mCountLeft += left;
    mCountRight += right;
    incrementBest(evidence, read);
  }

  /**
   * Increment the statistics.
   * @param evidence with probabilities for hypotheses.
   * @param bestHyp code for the best hypothesis supported by this piece of evidence
   */
  protected abstract void incrementBest(final EvidenceInterface evidence, int bestHyp);

  /**
   * @return the effective coverage for the model
   */
  public abstract int coverage();

  /**
   * @return the number of reads which did not equal the reference.
   */
  public final int nonRefCount() {
    return (int) MathUtils.round(mNonRefCount);
  }

  /**
   * @return the number of ambiguous reads.
   */
  public abstract int ambiguousCount();

  /**
   * @return the total of the error corrections.
   */
  public abstract double totalError();

  /**
   * The count of unmapped evidence seen at this locus.
   * @return unmapped evidence count
   */
  public int placedUnmappedCount() {
    return mCountUnmapped;
  }

  /**
   * The ratio of placed unmapped evidence to total placed evidence.
   * @return unmapped evidence ratio
   */
  public Double placedUnmappedRatio() {
    final int mappedCoverage = coverage();
    if (mappedCoverage > 0) {
      return placedUnmappedCount() / (double) (mappedCoverage + placedUnmappedCount());
    }
    return null;
  }

  /**
   * Replaces allele statistics on this instance with one remapping the
   * description into a new description according to specified description mappings.
   * @param newDescription the new description
   * @param mapping a mapping array. counts for old description <code>i</code> should
   * map to new description <code>mapping[i]</code>
   */
  public void remapAlleleStatistics(Description newDescription, int[] mapping) {
    setCounts(counts().remap(newDescription, mapping));
  }

  /**
   * Determine the ambiguity ratio at the current position
   * @return the ambiguity ratio, null if no counts.
   */
  public final Double ambiguityRatio() {
    final int count = coverage();
    if (count > 0) {
      return ambiguousCount() / (double) count;
    }
    return null;
  }

  /**
   * Check if the current position is over the coverage threshold.
   * @param params command line parameters.
   * @param sequenceName what sequence is the call in
   * @return true iff over coverage threshold.
   */
  public final boolean overCoverage(final VariantOutputOptions params, String sequenceName) {
    return coverage() > params.maxCoverageFilter().thresholdSingle(sequenceName);
  }

  /**
   * Determine if the current position is ambiguous.
   * @param params command line parameters.
   * @return true iff ambiguous.
   */
  public final boolean ambiguous(final VariantOutputOptions params) {
    //TODO use z-score to make more robust.
    final Double maxAmbiguity = params.maxAmbiguity();
    final Double ambiguityRatio = ambiguityRatio();
    //System.err.println("mAmbiguous=" + mAmbiguous + ", total=" + mCounts.totalCount() + ", ambiguityRatio=" + ambiguityRatio);
    return maxAmbiguity != null && ambiguityRatio != null && ambiguityRatio > maxAmbiguity;
  }

  /**
   * Decorate a <code>Variant</code> with information for a call.
   * @param variant the variant to decorate
   * @param templateName name of the template
   * @param params parameters controlling output options
   */
  public final void addFiltersToVariant(Variant variant, String templateName, VariantOutputOptions params) {
    //TODO shouldn't we set both of these if necessary?
    if (overCoverage(params, templateName)) {
      variant.addFilter(VariantFilter.COVERAGE);
    } else if (ambiguous(params)) {
      variant.addFilter(VariantFilter.AMBIGUITY);
    }
  }

  /**
   * Decorate a <code>VariantSample</code> with information for a call.
   *
   * @param sample variant sample object
   * @param model Model with details including posteriors.
   * @param params parameters controlling output options
   */
  public void addCountsToSample(VariantSample sample, ModelInterface<?> model, VariantOutputOptions params) {
    sample.setAmbiguityRatio(ambiguityRatio());
    sample.setCoverage(coverage());
    sample.setCoverageCorrection(totalError());
    if (sample.getName() != null) {
      final double unmatedProbability = unmatedProbability();
      if (sample.getPloidy() == Ploidy.HAPLOID) {
        final int alleleIndex = getAlleleIndex(sample.getName(), mCounts.getDescription());
        sample.setHoeffdingStrandBiasAllele1(mCounts.strandBias(alleleIndex));
        sample.setHoeffdingUnmatedBiasAllele1(mCounts.unmatedBias(alleleIndex, unmatedProbability));
        sample.setHoeffdingAlleleBalanceHom(mCounts.alleleBalanceHomozygous(alleleIndex, coverage()));
      } else if (sample.getPloidy() == Ploidy.DIPLOID) {
        getAlleleIndexes(sample.getName(), mCounts.getDescription(), mAlleleIndexes);
        sample.setHoeffdingStrandBiasAllele1(mCounts.strandBias(mAlleleIndexes[0]));
        sample.setHoeffdingStrandBiasAllele2(mCounts.strandBias(mAlleleIndexes[1]));
        sample.setHoeffdingUnmatedBiasAllele1(mCounts.unmatedBias(mAlleleIndexes[0], unmatedProbability));
        sample.setHoeffdingUnmatedBiasAllele2(mCounts.unmatedBias(mAlleleIndexes[1], unmatedProbability));
        if (mAlleleIndexes[0] == mAlleleIndexes[1]) {
          sample.setHoeffdingAlleleBalanceHom(mCounts.alleleBalanceHomozygous(mAlleleIndexes[0], coverage()));
        } else {
          sample.setHoeffdingAlleleBalanceHet(mCounts.alleleBalance(mAlleleIndexes[0], mAlleleIndexes[1]));
        }
      }
      sample.setHoeffdingReadPositionBias(readPositionBias());
    }
    sample.setStats(this);
    sample.setPlacedUnmappedRatio(placedUnmappedRatio());
  }

  protected abstract double unmatedProbability();

  protected Double readPositionBias() {
    final int trials = mCountLeft + mCountRight;
    final int observed = mCountLeft;
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  private static int getAlleleIndex(final String name, final Description des) {
    for (int j = 0; j < des.size(); ++j) {
      if (name.equals(des.name(j))) {
        return j;
      }
    }
    return -1;
  }

  private static void getAlleleIndexes(String names, Description des, int[] alleleIndexes) {
    final String[] split = StringUtils.split(names, VariantUtils.COLON);
    assert split.length == 2;
    Arrays.fill(alleleIndexes, -1);
    for (int i = 0; i < split.length; ++i) {
      for (int j = 0; j < des.size(); ++j) {
        if (split[i].equals(des.name(j))) {
          alleleIndexes[i] = j;
        }
      }
    }
//    if (alleleIndexes[0] == -1 || alleleIndexes[1] == -1) {
//      System.out.println(names);
//    }
  }

  /**
   * Underlying counts for each hypothesis.
   * @return counts
   */
  public T counts() {
    return mCounts;
  }

  public void setCounts(T counts) {
    mCounts = counts;
  }

  @Override
  public String toString() {
    // only used for debugging
    final StringBuilder sb = new StringBuilder();
    sb.append("Counts ");
    sb.append(StringUtils.LS);
    sb.append("coverage=").append(Integer.toString(coverage()));
    final String correction = String.format("%1$04.3f", totalError());
    sb.append(" correction=").append(correction);
    sb.append(StringUtils.LS);
    sb.append(mCounts);
    return sb.toString();
  }

  /**
   * Make a copy of this statistics object
   * @return copy
   */
  public Object copy() {
    try {
      return super.clone();
    } catch (CloneNotSupportedException e) {
      //actually it is
    }
    return null;
  }
}
