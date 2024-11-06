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

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Covariate;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.reader.Arm;
import com.rtg.util.Histogram;
import com.rtg.util.MathUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.VariantUtils;

/**
 * Implementation that uses calibration files.
 */
public class CalibratedMachineErrorParams extends AbstractMachineErrorParams {

  //any phred score lower than 2 theoretically means that all other options are more likely than the observed value.
  private static final int MIN_PHRED_SCORE = 2;

  private final PhredScaler mScaler;
  private final MachineType mMachineType;

  private final double[] mErrorDelDist;
  private final double[] mErrorInsDist;
  private final double[] mErrorMnpDist;
  private final double[] mErrorGapDist;
  private final double[] mErrorSmallGapDist;
  private final double[] mErrorOverlapDist;
  private final double[] mErrorOverlap2Dist;

  private final double mDelBaseRate;
  private final double mInsBaseRate;
  private final double mMnpBaseRate;
  private final double mDelEventRate;
  private final double mInsEventRate;
  private final double mMnpEventRate;

  /**
   * Constructor
   * @param mt machine type for read group
   * @param calibrator calibrator for all read groups
   * @param readGroupId read group to use for this error parameters.
   */
  public CalibratedMachineErrorParams(MachineType mt, Calibrator calibrator, String readGroupId) {
    mMachineType = mt;

    final Histogram nhHist = calibrator.getHistogram(Calibrator.NH_DIST, readGroupId);
    if (nhHist == null) {
      throw new NoTalkbackSlimException("No calibration data supplied for read group: " + readGroupId);
    } else if (nhHist.getLength() == 0) {
      throw new NoTalkbackSlimException("Calibration data indicates no coverage for the read group: " + readGroupId);
    }
    final Histogram delHist = calibrator.getHistogram(Calibrator.DEL_DIST, readGroupId);
    mErrorDelDist = delHist != null && delHist.getLength() > 1 ? delHist.toDistribution() : new double[] {0.0, 1.0};
    final Histogram insHist = calibrator.getHistogram(Calibrator.INS_DIST, readGroupId);
    mErrorInsDist = insHist != null && insHist.getLength() > 1 ? insHist.toDistribution() : new double[] {0.0, 1.0};
    final Histogram mnpHist = calibrator.getHistogram(Calibrator.MNP_DIST, readGroupId);
    mErrorMnpDist = mnpHist != null && mnpHist.getLength() > 1 ? mnpHist.toDistribution() : new double[] {0.0, 1.0};

    double[][] dists = null;
    if (isCG()) {
      if (mMachineType == MachineType.COMPLETE_GENOMICS_2) {
        if (calibrator.hasHistogram(Calibrator.CGOVER_DIST, readGroupId)) {
          final Histogram olapH = calibrator.getHistogram(Calibrator.CGOVER_DIST, readGroupId);
          dists = histogramsToCgV2Dists(olapH);
        }
      } else if (mMachineType == MachineType.COMPLETE_GENOMICS_2) {
        if (calibrator.hasHistogram(Calibrator.CGGAP_DIST, readGroupId)) {
          final Histogram gapH = calibrator.getHistogram(Calibrator.CGGAP_DIST, readGroupId);
          //overlap can be all 0s if we are parsing an old-style CG sam file, in which case use default overlap.
          final Histogram olapH;
          if (calibrator.hasHistogram(Calibrator.CGOVER_DIST, readGroupId)) {
            olapH = calibrator.getHistogram(Calibrator.CGOVER_DIST, readGroupId);
          } else {
            olapH = null;
          }
          dists = histogramsToCgV1Dists(gapH, olapH);
        }
      } else {
        throw new RuntimeException("CG calibration not supported for machine type: " + mMachineType);
      }
      if (dists == null) {
        Diagnostic.developerLog("Calibration file without CG gap/overlap information, using default priors");
      }
    }
    mErrorGapDist = dists != null && dists[0] != null ? dists[0] : MachineErrorParamsBuilder.CG_DEFAULT_GAP_DIST;
    mErrorSmallGapDist = dists != null && dists[1] != null ? dists[1] : MachineErrorParamsBuilder.CG_DEFAULT_SMALL_GAP_DIST;
    mErrorOverlapDist = dists != null && dists[2] != null ? dists[2] : MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP_DIST;
    mErrorOverlap2Dist = dists != null && dists[3] != null ? dists[3] : MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP2_DIST;

    final CalibrationStats sums = calibrator.getSums(CovariateEnum.READGROUP, readGroupId);
    //TODO Decide if Deletion or Insertion  should be included in potential error sites
    final long potentialErrorSites = sums.getEqual() + sums.getDifferent() + 1;
    mMnpBaseRate = (double) sums.getDifferent() / potentialErrorSites;
    mInsBaseRate = (double) sums.getInserted() / potentialErrorSites;
    mDelBaseRate = (double) sums.getDeleted() / potentialErrorSites;

    mMnpEventRate = mMnpBaseRate / GenomePriorParamsBuilder.averageLength(mErrorMnpDist);
    mInsEventRate = mInsBaseRate / GenomePriorParamsBuilder.averageLength(mErrorInsDist);
    mDelEventRate = mDelBaseRate / GenomePriorParamsBuilder.averageLength(mErrorDelDist);

    mScaler = getScaler(calibrator, readGroupId);
  }

  private static PhredScaler getScaler(Calibrator cal, String readGroup) {
    if (cal.getCovariateIndex(CovariateEnum.BASEQUALITY) == -1) {
      throw new IllegalArgumentException("Base quality covariate currently required");
    }
    final Calibrator.QuerySpec query = cal.initQuery();
    final int readGroupIndex = cal.getCovariateIndex(CovariateEnum.READGROUP);
    if (readGroupIndex >= 0) {
      query.setValue(CovariateEnum.READGROUP, cal.getCovariate(readGroupIndex).parse(readGroup));
    }
    if (cal.getCovariateIndex(CovariateEnum.ARM) != -1 && cal.getCovariateIndex(CovariateEnum.MACHINECYCLE) != -1) {
      final PhredScaler left = getArmPhredScaler(cal, query, Arm.LEFT);
      final PhredScaler right = getArmPhredScaler(cal, query, Arm.RIGHT);
      query.setValue(CovariateEnum.ARM, -1);
      return new ArmMachineCyclePhredScaler(left, right);
    } else if (cal.getCovariateIndex(CovariateEnum.MACHINECYCLE) != -1) {
      return new BaseQualityMachineCyclePhredScaler(cal, query);
    } else {
      return new BaseQualityPhredScaler(cal, query);
    }
  }

  private static PhredScaler getArmPhredScaler(Calibrator cal, Calibrator.QuerySpec query, Arm left) {
    query.setValue(CovariateEnum.ARM, left.ordinal());
    try {
      if (GlobalFlags.getBooleanValue(CoreGlobalFlags.QUALITY_CALIBRATION_COVARIATE_INTERSECTION)) {
        return new CovariateIntersectionCycleQualityScaler(cal, query);
      } else {
        return new BaseQualityMachineCyclePhredScaler(cal, query);
      }
    } finally {
      query.setValue(CovariateEnum.ARM, -1);
    }
  }

  static int countsToEmpiricalQuality(long mismatches, long totalReadPositions, double globalErrorRate) {
    final double error = (mismatches + globalErrorRate) / (totalReadPositions + 1.0);
    return Math.max((int) MathUtils.round(MathUtils.phred(error)), MIN_PHRED_SCORE);
  }

  /**
   * Produce gap/overlap distributions for CG version 1 data
   * @param gapHistogram histogram with index equal to gap length
   * @param overlapHistogram histogram with index equal to overlap length
   * @return the four distributions (gap, small gap, overlap, overlap 2)
   */
  static double[][] histogramsToCgV1Dists(Histogram gapHistogram, Histogram overlapHistogram) {
    //Partition and shunt histograms into expected array sizes and orientation
    //each read must have 1 and only 1 gap of between 4 and 8 (incl) long
    long totalReads = 0;
    for (int i = 4; i < gapHistogram.getLength(); ++i) {
      totalReads += gapHistogram.getValue(i);
    }

    //gap distribution expected with index 0 being probability of gap length 4 and index 4 being gap length 8;
    final double[] gapDist = new double[5];
    for (int i = 0; i < gapDist.length && 4 + i < gapHistogram.getLength(); ++i) {
      gapDist[i] = (double) gapHistogram.getValue(4 + i) / totalReads;
    }

    //each read must have 1 and only 1 gap of between 0 and 3 (incl) long, but we can't count those with 0 long gaps
    int totSmallGaps = 0;
    for (int i = 0; i < 4 && i < gapHistogram.getLength(); ++i) {
      totSmallGaps += gapHistogram.getValue(i);
    }
    if (totSmallGaps > totalReads) {
      throw new IllegalArgumentException("invalid CG gap distribution, small gaps (1-3) has larger count than big gaps(4-8)");
    }
    if (gapHistogram.getLength() > 9) {
      throw new IllegalArgumentException("invalid CG gap distribution");
    }

    //small gaps in slightly easier form of entry being probabilty of gap with length of its index
    double[] smallGapDist = new double[4];
    long totalSmall = 0;
    for (int i = 1; i < smallGapDist.length && i < gapHistogram.getLength(); ++i) {
      smallGapDist[i] = (double) gapHistogram.getValue(i) / totalReads;
      totalSmall += gapHistogram.getValue(i);
    }
    if (totalSmall == 0) {
      smallGapDist = MachineErrorParamsBuilder.CG_DEFAULT_SMALL_GAP_DIST;
    } else {
      //infer number of 0 long gaps from total reads and counted small gaps
      smallGapDist[0] = (double) (totalReads - totSmallGaps) / (double) totalReads;
    }

    final double[] overlapDist;
    if (overlapHistogram != null) {
      //each read must have 1 and only 1 overlap of between 0 and 4 (incl) long, again we can't keep track of 0 long whilst counting
      long totalOverlap = 0;
      for (int i = 1; i < overlapHistogram.getLength(); ++i) {
        totalOverlap += overlapHistogram.getValue(i);
      }
      if (totalOverlap > totalReads) {
        throw new IllegalArgumentException("invalid CG overlap distribution, more overlap than big gaps");
      }
      if (overlapHistogram.getLength() > 5) {
        throw new IllegalArgumentException("invalid CG overlap distribution");
      }

      //machine error params expects overlap distribution in wierd form with index 0 referring to probabilty of gap of 4 long and index 4 referring to probability of gap of 0 long
      overlapDist = new double[5];
      for (int i = 1; i < overlapHistogram.getLength(); ++i) {
        overlapDist[4 - i] = (double) overlapHistogram.getValue(i) / totalReads;
      }
      //again use total to infer number of 0 long overlaps
      overlapDist[4] = (double) (totalReads - totalOverlap) / (double) totalReads;
    } else {
      overlapDist = MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP_DIST;
    }
    return new double[][] {gapDist, smallGapDist, overlapDist, MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP2_DIST};
  }

  /**
   * Produce gap/overlap distributions for CG version 2 data
   * @param overlapHistogram histogram with index equal to overlap length
   * @return the four distributions (gap, small gap, overlap, overlap 2)
   */
  static double[][] histogramsToCgV2Dists(Histogram overlapHistogram) {
    //Shunt histogram into expected array size and orientation
    final double[] overlapDist;
    if (overlapHistogram != null) {
      //each read must have 1 and only 1 overlap of between 0 and 7 (incl) long, again we can't keep track of 0 long whilst counting
      long totalOverlap = overlapHistogram.getLength(); // For Laplace estimator
      for (int i = 1; i < overlapHistogram.getLength(); ++i) {
        totalOverlap += overlapHistogram.getValue(i);
      }
      if (overlapHistogram.getLength() > 8) {
        throw new IllegalArgumentException("invalid CG overlap distribution");
      }

      //machine error params expects overlap distribution in wierd form with index 0 referring to probabilty of gap of 7 long and index 7 referring to probability of gap of 0 long
      overlapDist = new double[8];
      for (int i = 0; i < overlapHistogram.getLength(); ++i) {
        overlapDist[7 - i] = (double) (overlapHistogram.getValue(i) + 1) / totalOverlap;
      }
    } else {
      overlapDist = MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP_DIST;
    }
    return new double[][] {MachineErrorParamsBuilder.CG_DEFAULT_GAP_DIST, MachineErrorParamsBuilder.CG_DEFAULT_SMALL_GAP_DIST, MachineErrorParamsBuilder.CG_DEFAULT_OVERLAP_DIST, overlapDist};
  }

  @Override
  public double errorDelBaseRate() {
    return mDelBaseRate;
  }

  @Override
  public double[] errorDelDistribution() {
    return mErrorDelDist;
  }

  @Override
  public double errorDelEventRate() {
    return mDelEventRate;
  }

  @Override
  public double errorInsBaseRate() {
    return mInsBaseRate;
  }

  @Override
  public double[] errorInsDistribution() {
    return mErrorInsDist;
  }

  @Override
  public double errorInsEventRate() {
    return mInsEventRate;
  }

  @Override
  public double[] errorMnpDistribution() {
    return mErrorMnpDist;
  }

  @Override
  public double errorMnpEventRate() {
    return mMnpEventRate;
  }

  @Override
  public double errorSnpRate() {
    return mMnpBaseRate;
  }

  @Override
  public double[] gapDistribution() {
    return mErrorGapDist;
  }


  @Override
  public double[] smallGapDistribution() {
    return mErrorSmallGapDist;
  }

  @Override
  public double[] overlapDistribution() {
    return mErrorOverlapDist;
  }

  @Override
  public double[] overlapDistribution2() {
    return mErrorOverlap2Dist;
  }

  @Override
  public int getScaledPhred(byte quality, int readPos, Arm arm) {
    return mScaler.getScaledPhred(quality, readPos, arm);
  }

  @Override
  public final boolean isCG() {
    return mMachineType.isCG();
  }

  @Override
  public MachineType machineType() {
    return mMachineType;
  }


  @Override
  protected int[] qualityCurve() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  private static final String DUMP_MERGED = "dump-merged";
  private static final String DUMP_PARAMS = "dump-params";

  /**
   * Test from command line
   * @param args command line arguments
   * @throws IOException if the calibration files cannot be read
   */
  public static void main(String... args) throws IOException {
    final CFlags flags = new CFlags();
    flags.registerOptional(DUMP_MERGED, "if set, output merged calibration file");
    flags.registerOptional(DUMP_PARAMS, "if set, output machine error params");
    flags.registerRequired(File.class, CommonFlags.FILE, "SAM file containing unmapped reads to attempt to rescue").setMaxCount(Integer.MAX_VALUE);
    flags.setFlags(args);
    Diagnostic.setLogStream();
    Covariate[] covariates = null;
    for (final Object o : flags.getAnonymousValues(0)) {
      covariates = Calibrator.getCovariateSet((File) o);
    }
    if (covariates == null) {
      throw new IllegalStateException("no calibration covariates found");
    }
    final Calibrator c = new Calibrator(covariates, null);
    for (final Object o : flags.getAnonymousValues(0)) {
      c.accumulate((File) o);
    }
    if (flags.isSet(DUMP_MERGED)) {
      c.writeToStream(System.out);
    }
    if (flags.isSet(DUMP_PARAMS)) {
      final Covariate rgc = c.getCovariate(c.getCovariateIndex(CovariateEnum.READGROUP));
      for (int rgid = 0; rgid < rgc.size(); ++rgid) {
        System.out.println("##############################################");
        final String rgname = rgc.valueString(rgid);
        final CalibratedMachineErrorParams cme = new CalibratedMachineErrorParams(MachineType.ILLUMINA_PE, c, rgname);
        System.out.println("# Machine errors for read group: " + rgname);
        System.out.println(VariantUtils.toMachineErrorProperties(cme));
      }
    }
  }
}
