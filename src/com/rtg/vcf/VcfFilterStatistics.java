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

package com.rtg.vcf;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.launcher.Statistics;

/**
 * Encapsulate the counts of filtered VCF records
 */
class VcfFilterStatistics implements Statistics {

  enum Stat {
    SAME_AS_REF_FILTERED_COUNT
    , ALL_SAME_AS_REF_FILTERED_COUNT
    , AMBIGOUS_FILTERED_COUNT
    , READ_DEPTH_FILTERED_COUNT
    , FAILED_KEEP_COUNT
    , NOT_SNP_COUNT
    , GENOTYPE_QUALITY_POSTERIOR_FILTERED_COUNT
    , QUALITY_FILTERED_COUNT
    , ALLELE_BALANCE_FILTERED_COUNT
    , DENSITY_WINDOW_COUNT
    , EXCLUDE_BED_COUNT
    , INCLUDE_BED_COUNT
    , WRITTEN_COUNT
    , TOTAL_COUNT
    , AVR_SCORE_FILTERED_COUNT
    , OVERLAP_COUNT
    , SNP_COUNT
    , DENOVO_SCORE
    , COMBINED_READ_DEPTH_FILTERED_COUNT
  }

  private int[] mValues = new int[Stat.values().length];
  private int[] mFilteredCount;
  private int[] mInfoCount;
  final Map<String, Integer> mFilterTags = new TreeMap<>(); // valid filter tags - input file specific
  final Map<String, Integer> mInfoTags = new TreeMap<>(); // valid info tags - input file specific
  private boolean mPosteriorFiltering;

  @Override
  public void printStatistics(final OutputStream stream) {
    if (stream != null) {
      final PrintStream output = new PrintStream(stream);
      output.println();
      output.println("Total records : " + get(Stat.TOTAL_COUNT));

      for (final String tag : mFilterTags.keySet()) {
        printCount(output, tag, mFilteredCount[mFilterTags.get(tag)]);
      }
      for (final String tag : mInfoTags.keySet()) {
        printCount(output, tag, mInfoCount[mInfoTags.get(tag)]);
      }
      printCount(output, "quality", get(Stat.QUALITY_FILTERED_COUNT));
      if (isPosteriorFiltering()) {
        printCount(output, "posterior", get(Stat.GENOTYPE_QUALITY_POSTERIOR_FILTERED_COUNT));
      } else {
        printCount(output, "genotype quality", get(Stat.GENOTYPE_QUALITY_POSTERIOR_FILTERED_COUNT));
      }
      printCount(output, "AVR score", get(Stat.AVR_SCORE_FILTERED_COUNT));
      printCount(output, "sample read depth", get(Stat.READ_DEPTH_FILTERED_COUNT));
      printCount(output, "combined read depth", get(Stat.COMBINED_READ_DEPTH_FILTERED_COUNT));
      printCount(output, "ambiguity ratio", get(Stat.AMBIGOUS_FILTERED_COUNT));
      printCount(output, "allele balance", get(Stat.ALLELE_BALANCE_FILTERED_COUNT));
      printCount(output, "same as reference", get(Stat.SAME_AS_REF_FILTERED_COUNT));
      printCount(output, "all samples same as reference", get(Stat.ALL_SAME_AS_REF_FILTERED_COUNT));
      printCount(output, "not a SNP", get(Stat.NOT_SNP_COUNT));
      printCount(output, "simple SNP", get(Stat.SNP_COUNT));
      printCount(output, "not in keep set", get(Stat.FAILED_KEEP_COUNT));
      printCount(output, "overlap with previous", get(Stat.OVERLAP_COUNT));
      printCount(output, "density window", get(Stat.DENSITY_WINDOW_COUNT));
      printCount(output, "include file", get(Stat.INCLUDE_BED_COUNT));
      printCount(output, "exclude file", get(Stat.EXCLUDE_BED_COUNT));
      printCount(output, "de novo score", get(Stat.DENOVO_SCORE));
      //output.println("_Filtered due to other : " + m_Other_Filtered_COUNT);
      output.println("Remaining records : " + get(Stat.WRITTEN_COUNT));

      //assert mTotalCount == (mComplexFilteredCount + mDensityFilteredCount + mGenotypeQualityPosteriorFilteredCount + mReadDepthFilteredCount + mWrittenCount + mAmbigousFilteredCount + mSameAsRefFilteredCount + mNotSnpCount + mAlleleBalanceFilteredCount + mQualityFilteredCount + mComplexRegionFilteredCount + mOtherFilteredCount);
    }
  }

  @Override
  public void generateReport() {
    // Unimplemented
  }


  private void printCount(PrintStream output, String tag, int count) {
    if (count != 0) {
      output.println("Filtered due to " + tag + " : " + count);
    }
  }
  int get(Stat s) {
    return mValues[s.ordinal()];
  }
  public void increment(Stat s) {
    mValues[s.ordinal()]++;
  }
  void incrementFilterTag(String tag) {
    mFilteredCount[mFilterTags.get(tag)]++;
  }
  void incrementInfoTag(String tag) {
    mInfoCount[mInfoTags.get(tag)]++;
  }

  boolean isPosteriorFiltering() {
    return mPosteriorFiltering;
  }

  void setPosteriorFiltering(boolean posteriorFiltering) {
    mPosteriorFiltering = posteriorFiltering;
  }

  void setFilterTags(Map<String, Integer> filterTags) {
    mFilterTags.putAll(filterTags);
    mFilteredCount = new int[mFilterTags.size()];
  }

  void setInfoTags(Map<String, Integer> infoTags) {
    mInfoTags.putAll(infoTags);
    mInfoCount = new int[mInfoTags.size()];
  }
}
