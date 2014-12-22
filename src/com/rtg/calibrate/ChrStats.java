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
package com.rtg.calibrate;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.Set;

import com.rtg.reader.SequencesReader;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Checks whether chromosome coverage levels are consistent with expected levels for the specified sex.
 */
public class ChrStats {

  private static final long MIN_SPECIFIED_SEQUENCES = 5;
  static final double DEFAULT_MIN_SEX_DEVIATIONS = 5;
  static final double DEFAULT_MIN_DEVIATIONS = 10;
  private static final String INCONSISTENT = " ** INCONSISTENT **";
  private static final String SEX_INCONSISTENT = " ** SEX INCONSISTENT **";

  /** Hold the result of a chromosome statistics check. */
  public static final class ChrStatsResult {
    private final Sex mClaimedSex;
    private boolean mIsSexConsistent;
    private final int mStrangeCount;
    private final int mTotalCount;
    private final double mMeanDipCoverage;
    private Sex mObservedSex;

    private ChrStatsResult(final Sex claimedSex, final boolean isSexConsistent, final int strangeCount, final int totalCount, final double meanDipCoverage) {
      mClaimedSex = claimedSex;
      mObservedSex = null;
      mIsSexConsistent = isSexConsistent;
      mStrangeCount = strangeCount;
      mTotalCount = totalCount;
      mMeanDipCoverage = meanDipCoverage;
    }

    public Sex getClaimedSex() {
      return mClaimedSex;
    }

    public Sex getObservedSex() {
      return mObservedSex;
    }
  }


  static final NumberFormat NF = NumberFormat.getInstance();
  static {
    NF.setMaximumFractionDigits(2);
    NF.setMinimumFractionDigits(2);
  }

  private final SequencesReader mGenomeReader;
  private final double mMinSexDeviations;
  private final double mMinDeviations;

  /**
   * Construct a new coverage checker reporting at a given number of deviations.
   * @param genomeReader genome related information
   * @param sexDeviations deviations needed before a potential sex problem is reported
   * @param deviations deviations needed to count as bad
   */
  public ChrStats(final SequencesReader genomeReader, final double sexDeviations, double deviations) {
    mGenomeReader = genomeReader;
    mMinSexDeviations = sexDeviations;
    mMinDeviations = deviations;
  }

  /**
   * Construct a new coverage checker.
   * @param genomeParams genome related information
   */
  public ChrStats(final SequencesReader genomeParams) {
    this(genomeParams, DEFAULT_MIN_SEX_DEVIATIONS, DEFAULT_MIN_DEVIATIONS);
  }

  /**
   * Compare the coverage of each "specified" sequence in the reference against the average of the
   * coverage for all the other "specified" sequences in the reference.  Significant departures for
   * are reported.
   * @param calibrator calibrator containing coverage information to use
   * @param sample names of sample
   * @param sex claimed sex for the sample
   * @param log should the results of the calculation be logged
   * @return true iff the sample has a coverage anomaly when treated as the specified sex
   */
  ChrStatsResult coverageCheck(final CalibratedPerSequenceExpectedCoverage calibrator, final String sample, final Sex sex, final boolean log) throws IOException {
    Diagnostic.userLog("Starting check of " + sample + " for sex=" + sex + " with autosome-z-threshold=" + mMinDeviations + " and sex-z-threshold=" + mMinSexDeviations);
    // If the fallback fires below, then the test will not end up happening because the number
    // of specified sequences will be 0.  But this way prevents an exception from occurring if
    // there is no reference.txt file present in the template SDF.
    final ReferenceGenome referenceGenome = new ReferenceGenome(mGenomeReader, sex, ReferenceGenome.DefaultFallback.DIPLOID);
    int numSpecifiedSequences = 0;
    double sumCoverage = 0;
    double sumCoverageSquares = 0;
    // Compute total sum
    for (final ReferenceSequence seq : referenceGenome.sequences()) {
      if (seq.isSpecified() && seq.ploidy() == Ploidy.DIPLOID) {
        final double coverage = calibrator.expectedCoverage(seq.name(), sample);
        sumCoverage += coverage;
        sumCoverageSquares += coverage * coverage;
        numSpecifiedSequences++;
      }
    }
    if (numSpecifiedSequences < MIN_SPECIFIED_SEQUENCES) {
      return new ChrStatsResult(sex, true, 0, numSpecifiedSequences, 0); // Not enough specified sequences to reliably make these checks
    }
    final double meanDip = sumCoverage / numSpecifiedSequences;
    final double varDip = sumCoverageSquares / numSpecifiedSequences - meanDip * meanDip;
    if (varDip < 0) {
      Diagnostic.userLog("Arithmetic underflow " + sample + " mean=" + meanDip + " var=" + varDip);
      return new ChrStatsResult(sex, true, 0, numSpecifiedSequences, meanDip);
    }
    // Iterate again over sequences doing leave one out for diploid chromosomes or
    // corrected comparison with all diploid coverage in the case of haploid chromosomes
    boolean consistent = true;
    final double stddevDip = Math.sqrt(varDip);
    final long leaveOneOut = numSpecifiedSequences - 1;
    int strangeCount = 0;
    int totalCount = 0;
    for (final ReferenceSequence seq : referenceGenome.sequences()) {
      if (seq.isSpecified()) {
        totalCount++;
        final double coverage = calibrator.expectedCoverage(seq.name(), sample);
        if (seq.ploidy() == Ploidy.DIPLOID) {
          final double mean = (sumCoverage - coverage) / leaveOneOut;
          final double var = (sumCoverageSquares - coverage * coverage) / leaveOneOut - mean * mean;
          if (var < 0) {
            Diagnostic.userLog("Arithmetic underflow " + sample + " mean=" + mean + " var=" + var + " seq=" + seq.name());
            continue;
          }
          final double stddev = Math.sqrt(var);
          final double z = (coverage - mean) / stddev;
          String status = "";
          final double abs = Math.abs(z);
          if (abs > mMinSexDeviations) {
            status = SEX_INCONSISTENT;
            consistent = false;
          }
          if (abs > mMinDeviations) {
            status = INCONSISTENT;
            strangeCount++;
          }
          if (log) {
            Diagnostic.userLog(sample + " " + seq.name() + " observed-coverage=" + NF.format(coverage) + " z=" + NF.format(z) + "  (diploid coverage mean=" + NF.format(mean) + " stddev=" + NF.format(stddev) + ")" + status);
          }
        } else if (seq.ploidy() == Ploidy.HAPLOID) {
          final double z = (2 * coverage - meanDip) / stddevDip;
          String status = "";
          final double abs = Math.abs(z);
          if (abs > mMinSexDeviations) {
            status = SEX_INCONSISTENT;
            consistent = false;
          }
          if (abs > mMinDeviations) {
            status = INCONSISTENT;
            strangeCount++;
          }
          if (log) {
            Diagnostic.userLog(sample + " " + seq.name() + " observed-coverage=" + NF.format(coverage) + " z=" + NF.format(z) + "  (diploid coverage mean=" + NF.format(meanDip) + " stddev=" + NF.format(stddevDip) + ")" + status);
          }
        }
      }
    }
    return new ChrStatsResult(sex, consistent, strangeCount, totalCount, meanDip);
  }

  /**
   * Run a chromosome statistics check for a sample.
   * @param calibrator calibrator containing coverage information to use
   * @param sample names of sample
   * @param sex claimed sex for the sample
   * @return results
   * @throws IOException if an I/O error occurs.
   */
  public ChrStatsResult runCheckAndReport(final CalibratedPerSequenceExpectedCoverage calibrator, final String sample, final Sex sex) throws IOException {
    if (calibrator.containsSample(sample)) {
      // First run with the claimed sex including reporting of anamolous results
      final ChrStatsResult res = coverageCheck(calibrator, sample, sex, true);
      if (sex == Sex.EITHER) {
        // Try and guess the sex when none was specified, but suppress reporting of anamolous results
        final ChrStatsResult male = coverageCheck(calibrator, sample, Sex.MALE, false);
        final ChrStatsResult female = coverageCheck(calibrator, sample, Sex.FEMALE, false);
        if (male.mIsSexConsistent ^ female.mIsSexConsistent) {
          // Note male == true indicates pass male sex test, hence female and vice versa
          res.mObservedSex = male.mIsSexConsistent ? Sex.MALE : Sex.FEMALE;
          res.mIsSexConsistent = false;
        }
      } else if (!res.mIsSexConsistent) {
        // Sample fails with expected sex, swap sex and try again
        final Sex altSex = sex == Sex.MALE ? Sex.FEMALE : Sex.MALE;
        if (coverageCheck(calibrator, sample, altSex, false).mIsSexConsistent) {
          res.mObservedSex = altSex;
        }
        // else: could be a trisomy, XXY, etc., up to user to inspect logged results
      }
      return res;
    } else {
      Diagnostic.userLog(sample + " has no calibration information, no coverage checking done");
      return null;
    }
  }

  /**
   * String representation of sex.
   * @param sex the sex
   * @return string representation
   */
  public static String sexString(final Sex sex) {
    return sex == null ? "" : sex.toString();
  }

  private static TextTable initTable() {
    final TextTable table = new TextTable();
    table.setAlignment(TextTable.Align.LEFT, TextTable.Align.LEFT, TextTable.Align.LEFT, TextTable.Align.LEFT, TextTable.Align.RIGHT, TextTable.Align.LEFT);
    table.addRow("sample", "mapped-as", "consistent", "possible", "coverage", "");
    table.addSeparator();
    return table;
  }

  private static void addRow(final TextTable table, final String sample, final Sex sex, final ChrStatsResult res) {
    if (res == null) {
      table.addRow(sample, sexString(sex), "", "", "", "(no calibration information)");
    } else {
      final String countMessage = " (" + res.mStrangeCount + " of " + res.mTotalCount + " sequences have unexpected coverage level)";
      table.addRow(sample, sexString(sex), Boolean.toString(res.mIsSexConsistent), res.mIsSexConsistent ? "" : sexString(res.mObservedSex), Utils.realFormat(res.mMeanDipCoverage, 2), res.mStrangeCount > 0 ? countMessage : "");
    }
  }

  /**
   * Compare the coverage of each "specified" sequence in the reference against the average of the
   * coverage for all the other "specified" sequences in the reference.  Significant departures for
   * are reported.
   * @param calibrator calibrator containing coverage information to use
   * @param sample names of sample
   * @param sex claimed sex for the sample
   * @throws java.io.IOException if an I/O error occurs
   */
  public void chrStatsCheckAndReport(final CalibratedPerSequenceExpectedCoverage calibrator, final String sample, final Sex sex) throws IOException {
    final TextTable table = initTable();
    addRow(table, sample, sex, runCheckAndReport(calibrator, sample, sex));
    Diagnostic.warning(table.toString());
  }

  /**
   * Perform a basic leave one out coverage check for multiple samples.
   * @param calibrator calibrator containing coverage information to use
   * @param samples names of samples
   * @param pedigree pedigree information
   */
  void chrStatsCheck(final CalibratedPerSequenceExpectedCoverage calibrator, final Set<String> samples, final GenomeRelationships pedigree) throws IOException {
    final TextTable table = initTable();
    for (final String sample : samples) {
      final Sex sex = pedigree.getSex(sample);
      addRow(table, sample, sex, runCheckAndReport(calibrator, sample, sex));
    }
    System.out.print(table.toString());
  }
}
