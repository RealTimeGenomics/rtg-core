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
package com.rtg.variant.cnv;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;

import com.rtg.launcher.OutputParams;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.variant.cnv.region.EmptyRegion;
import com.rtg.variant.cnv.region.Region;
import com.rtg.variant.cnv.region.RegionUtils;

/**
 * Use two arrays of block counts and try a piece-wise constant fit to the data.
 *
 * Write out three files.
 * <code>directory/cnv.bed</code> with one line for each block.
 * Each line has one of the formats:
 * <ul>
 * <li><code>seq-id start end cnv mean standard-deviation</code>
 * <li><code>seq-id start end nregion</code>
 * <li><code>seq-id start end gdelregion[&lt;|&gt;]</code>
 * <li><code>seq-id start end model[0-5] mean</code>
 * </ul>
 * start is zero-based, end is exclusive, both in nucleotide space
 *  <code>directory/cnv.ratio</code>
 * One ratio number per line.
 */
public final class CnvRatio {

  private static final byte[] VERSION = "#Version ".getBytes();
  private static final byte[] CNV_VERSION = ", cnv v0.2".getBytes();
  private static final byte[] CL = "#CL".getBytes();
  private static final byte[] INTERNAL = "Internal".getBytes();
  private static final byte[] HEADERS_BED = (
      "#Seq\tstart\tend\tcnv\tmean\tstandard-deviation" + StringUtils.LS
      + "#Seq\tstart\tend\tnregion" + StringUtils.LS
      + "#Seq\tstart\tend\tgdelregion" + StringUtils.LS
      + "#Seq\tstart\tend\tmodel[0-5] mean" + StringUtils.LS).getBytes();
  private static final byte[] HEADERS_RATIO = ("#Seq\tposition\tratio" + StringUtils.LS).getBytes();
  private static final byte[] TAB = "\t".getBytes();
  private static final byte[] LS = StringUtils.LS.getBytes();
  private static final byte[] NREGION = "nregion".getBytes();
  private static final byte[] MODEL = "model".getBytes();
  private static final byte[] DIVDREGION = "gdelregion<".getBytes();
  private static final byte[] MULDREGION = "gdelregion>".getBytes();

  enum State {
    IN, OUT
  }

  private int mCnt1Total;

  private int mCnt2Total;

  private int mBlocks;

  private int mSeqStart;

  private State mState;

  private LSMOutput mSimpleOut;

  private OutputStream mRatio;

  private OutputStream mCnvOut;

  private LeastSquaresModel mLSM;

  private Region mNRegionExcludes = null;

  private Region mDivDelRegionExcludes = null;

  private Region mMulDelRegionExcludes = null;

  private String mSeqId;

  private final double mConstant;

  private final Map<String, Region> mNRegions;

  private final Map<Integer, String> mTemplateNameMap;

  private final OutputParams mOutputParams;

  private final int mBucketSize;

  private final double mMulFactor;

  private final double mDivFactor;
  private final double mMeanReadLength;
  private final boolean mExtraPenaltyOff;

  /**
   * CnvRatio constructor
   * @param constant reluctance to break
   * @param nRegions N-Region exclusion map
   * @param nameMap map of index to template names
   * @param params CNV product parameters
   * @param meanReadLength mean read length
   * @param penaltyOn flag if penalty is switched on
   */
  public CnvRatio(final double constant, final Map<String, Region> nRegions, final Map<Integer, String> nameMap,
      CnvProductParams params, final double meanReadLength, final boolean penaltyOn) {
    mConstant = constant;
    mNRegions = nRegions;
    mTemplateNameMap = nameMap;
    mOutputParams = params.outputParams();
    mBucketSize = params.bucketSize();
    mMulFactor = params.multiplicationFactor();
    mDivFactor = params.divisionFactor();
    mMeanReadLength = meanReadLength;
    mExtraPenaltyOff = penaltyOn;
  }

  /**
   * Run the ratio calculation over given data
   * @param refLine chunk counts for the base line
   * @param testLine chunk counts for the target
   * @throws IOException IO Exceptions happen sometimes
   */
  public void exec(final int[][] refLine, final int[][] testLine) throws IOException {
    mCnt1Total = 0;
    mCnt2Total = 0;
    assert refLine.length == testLine.length;
    try {
      mCnvOut = mOutputParams.outStream("cnv.bed");
      mRatio = mOutputParams.outStream("cnv.ratio");
      writeVersionCommandLine(mCnvOut);
      writeVersionCommandLine(mRatio);
      mCnvOut.write(HEADERS_BED);
      mRatio.write(HEADERS_RATIO);
      mSimpleOut = new SimpleOutput(mCnvOut, mBucketSize);
      for (int i = 0; i < refLine.length; i++) {
        //If not using all files in this analysis this can be true
        if (refLine[i] != null && testLine[i] != null) {
          assert refLine[i].length == testLine[i].length;
          startSeq(mTemplateNameMap.get(i), refLine[i]);
          for (int j = 0; j < refLine[i].length; j++) {
            add(refLine[i][j], testLine[i][j]);
          }
          endSeq(refLine[i].length);
        }
      }
      for (int i = 0; i < refLine.length; i++) {
        if (refLine[i] != null && testLine[i] != null) {
          outputModel(refLine[i].length, mTemplateNameMap.get(i));
        }
      }
    } finally {
      closeOutputs();
    }
  }

  void writeVersionCommandLine(OutputStream out) throws IOException {
    out.write(VERSION);
    out.write(Environment.getVersion().getBytes());
    out.write(CNV_VERSION);
    out.write(LS);
    out.write(CL);
    out.write(TAB);
    if (CommandLine.getCommandLine() != null) {
      out.write(CommandLine.getCommandLine().getBytes());
    } else {
      out.write(INTERNAL);
    }
    out.write(LS);
  }

  void startSeq(final String seqId, final int[] refLine) {
    mBlocks = 0;
    mSeqStart = -1;
    mState = State.OUT;
    mSeqId = seqId;
    mNRegionExcludes = mNRegions.containsKey(seqId) ? mNRegions.get(seqId) : EmptyRegion.EMPTY_REGION;
    mDivDelRegionExcludes = RegionUtils.findGermlineDeletesUnderMean(refLine, mDivFactor, mNRegionExcludes, mExtraPenaltyOff);
    mMulDelRegionExcludes = RegionUtils.findGermlineDeletesOverMean(refLine, mMulFactor, mNRegionExcludes, mExtraPenaltyOff);
    mLSM = new LeastSquaresModel(mSeqId, 2, mMeanReadLength / mBucketSize, mConstant, mExtraPenaltyOff, mSimpleOut);
  }

  void add(final int cnt1, final int cnt2) throws IOException {
    final double ratio = 2.0 * (cnt2 + 1.0) / (cnt1 + 1.0);
    final boolean okRegion = !mNRegionExcludes.isInRegion(mBlocks + 1) && !mDivDelRegionExcludes.isInRegion(mBlocks + 1) && !mMulDelRegionExcludes.isInRegion(mBlocks + 1);
    if (mState == State.IN) {
      if (!okRegion) {
        if (mLSM != null) {
          mLSM.scan(mSeqStart, mBlocks + 1);
        }
        mSeqStart = -1;
        mState = State.OUT;
      }
    } else if (okRegion) {
      mSeqStart = mBlocks + 1;
      mState = State.IN;
    }
    mLSM.add(ratio);
    if (okRegion) {
      mRatio.write(mSeqId.getBytes());
      mRatio.write(TAB);
      mRatio.write(Integer.toString(ntPos(mBlocks + 1)).getBytes());
      mRatio.write(TAB);
      mRatio.write(Utils.realFormat(ratio, 3).getBytes());
      mRatio.write(LS);
      mCnt1Total += cnt1;
      mCnt2Total += cnt2;
    }
    mBlocks++;
  }

  void endSeq(int seqSize) throws IOException {
    if (mLSM != null && mState == State.IN) {
      mLSM.scan(mSeqStart, mBlocks + 1);
    }
    outputRegions(seqSize, mNRegionExcludes, NREGION);
    outputRegions(seqSize, mDivDelRegionExcludes, DIVDREGION);
    outputRegions(seqSize, mMulDelRegionExcludes, MULDREGION);
  }

  void outputModel(int seqSize, String seqId) throws IOException {
    final double ratio = 2.0 * (mCnt2Total + 1.0) / (mCnt1Total + 1.0);
    final double gamma = 0.75;
    for (int i = 0; i <= 5; i++) {
      final double level = ratio * gamma + (1.0 - gamma) * i;
      final String ls = Utils.realFormat(level, 3);
      mCnvOut.write(seqId.getBytes());
      mCnvOut.write(TAB);
      mCnvOut.write('0');
      mCnvOut.write(TAB);
      mCnvOut.write(Integer.toString(ntPos(seqSize + 1)).getBytes());
      mCnvOut.write(TAB);
      mCnvOut.write(MODEL);
      mCnvOut.write(Integer.toString(i).getBytes());
      mCnvOut.write(TAB);
      mCnvOut.write(ls.getBytes());
      mCnvOut.write(LS);
    }

  }

  private void outputRegions(int seqSize, Region exclusions, byte[] label) throws IOException {
    int start = 0;
    boolean inRegion = false;
    for (int i = 1; i <= seqSize; i++) {
      final boolean exclude = exclusions.isInRegion(i);
      if (!inRegion && exclude) {
        start = i;
        inRegion = true;
      } else if (inRegion && !exclude) {
        writeRegion(start, i, label);
        inRegion = false;
      }
    }
    if (inRegion) {
      writeRegion(start, seqSize + 1, label);
    }
  }

  private void writeRegion(int start, int end, byte[] label) throws IOException {
    mCnvOut.write(mSeqId.getBytes());
    mCnvOut.write(TAB);
    mCnvOut.write(Integer.toString(ntPos(start)).getBytes());
    mCnvOut.write(TAB);
    mCnvOut.write(Integer.toString(ntPos(end)).getBytes());
    mCnvOut.write(TAB);
    mCnvOut.write(label);
    mCnvOut.write(LS);
  }

  private int ntPos(int blockPos) {
    return (blockPos - 1) * mBucketSize;
  }

  @SuppressWarnings("try")
  void closeOutputs() throws IOException {
    try (OutputStream ignored = mRatio; OutputStream ignored2 = mCnvOut) {
      mRatio = null;
      mCnvOut = null;
    }
  }
}
