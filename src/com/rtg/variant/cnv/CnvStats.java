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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;


/**
 */
public class CnvStats {

  int mNumFalsePositive;

  int mNumFalseNegative;

  int mNumCorrect;

  boolean mCollapse;

  private Writer mCorrectOut;
  private Writer mFalsePositiveOut;
  private Writer mFalseNegativeOut;

  static final int SEQUENCE_COLUMN = 0;
  static final int POSITION_START_COLUMN = 1;
  static final int POSITION_END_COLUMN = 2;
  static final int TYPE_COLUMN = 3;
  static final int MEAN_COLUMN = 4;

  static final String TB = "\t";
  private int mThreshold = 100;

  static class CnvRegionLine {
    String mSequence;
    int mPositionStart;
    int mPositionEnd;
    double mMean;
    String mType;

    CnvRegionLine(final String line) {
      final String[] parts = line.split(TB);
      mSequence = parts[SEQUENCE_COLUMN];
      mPositionStart = Integer.parseInt(parts[POSITION_START_COLUMN]);
      mPositionEnd = Integer.parseInt(parts[POSITION_END_COLUMN]);
      mType = parts[TYPE_COLUMN];
      if ("cnv".equals(mType)) {
        mMean = Double.parseDouble(parts[MEAN_COLUMN]);
      }
    }

    @Override
    public String toString() {
      return ""
          + mSequence + TB
          + mPositionStart + TB
          + mPositionEnd + TB
          + mType + TB
          + Utils.realFormat(mMean, 4)
        ;
    }
  }

  public void setThreshold(final int threshold) {
    this.mThreshold = threshold;
  }

  public void setCollapse(final boolean val) {
    mCollapse = val;
  }

  /**
   * set output files
   * @param outdir directory for output files
   */
  public void setOutput(final File outdir) {
    final File correctFile = new File(outdir, "correct");
    final File falsePositiveFile = new File(outdir, "false-positive");
    final File falseNegativeFile = new File(outdir, "false-negative");

    try {
      mCorrectOut = new OutputStreamWriter(FileUtils.createOutputStream(correctFile, false, false));
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_CREATED, correctFile.getPath());
    }
    try {
      mFalsePositiveOut = new OutputStreamWriter(FileUtils.createOutputStream(falsePositiveFile, false, false));
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_CREATED, falsePositiveFile.getPath());
    }
    try {
      mFalseNegativeOut = new OutputStreamWriter(FileUtils.createOutputStream(falseNegativeFile, false, false));
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_CREATED, falseNegativeFile.getPath());
    }
  }

  /**
   * close output files
   * @throws IOException If an IO error occurs
   */
  public void closeOutput() throws IOException {
    if (mCorrectOut != null) {
      mCorrectOut.close();
    }
    if (mFalsePositiveOut != null) {
      mFalsePositiveOut.close();
    }
    if (mFalseNegativeOut != null) {
      mFalseNegativeOut.close();
    }
  }

  private static String getLine(final BufferedReader r) throws IOException {
    String line = r.readLine();
    while (line != null && line.startsWith("#")) {
      line = r.readLine();
    }
    return line;
  }

  private static CnvRegionLine readRegionLine(final BufferedReader reader) throws IOException {
    String line;
    while ((line = getLine(reader)) != null) {
      //System.err.println(line);
      final CnvRegionLine region = new CnvRegionLine(line);
      if ("cnv".equals(region.mType)) {
        return region;
      }
    }
    return null;
  }

  static class RegionPosition {
    CnvRegionLine mRegion;
    private int mPos;
    private final List<CnvRegionLine> mList;
    private boolean mIsStart;

    RegionPosition(final List<CnvRegionLine> list) {
      mList = list;
      mIsStart = false;
      mPos = 0;
      internalNext();
    }

    int getPos() {
      return mIsStart ? mRegion.mPositionStart : mRegion.mPositionEnd;
    }
    boolean valid() {
      return mRegion != null;
    }
    String getSequence() {
      return mRegion.mSequence;
    }
    private void internalNext() {
      if (!mIsStart) {
        if (mPos < mList.size()) {
          mRegion = mList.get(mPos++);
        } else {
          mRegion = null;
        }
      }
      mIsStart = !mIsStart;
    }
    void next() {
      final int pos = getPos();
      while (valid() && pos == getPos()) {
        internalNext();
      }
    }
  }

  ArrayList<CnvRegionLine> readCnvRegionLines(final BufferedReader reader) throws IOException {
    final ArrayList<CnvRegionLine> ret = new ArrayList<>();
    CnvRegionLine current;
    while ((current = readRegionLine(reader)) != null) {
      ret.add(current);
    }
    return ret;
  }

  /**
   * Removes points that lie in the center of a group of regions with the same copy number.
   */
  void collapseCnvRegionLines(final ArrayList<CnvRegionLine> list) {
    int i = 0;
    while (i < list.size() - 1) {
      final CnvRegionLine current = list.get(i);
      final CnvRegionLine next = list.get(i + 1);
      if (Double.doubleToLongBits(current.mMean) == Double.doubleToLongBits(next.mMean) && current.mSequence.equals(next.mSequence)) {
        current.mPositionEnd = next.mPositionEnd;
        list.remove(i + 1);
      } else {
        i++;
      }
    }
  }

  void getStats(final Reader generatedReader, final Reader detectedReader, final Writer out) throws IOException {
    final ArrayList<CnvRegionLine> detected;
    try (BufferedReader bufDet = new BufferedReader(detectedReader)) {
      detected = readCnvRegionLines(bufDet);
    }
    final ArrayList<CnvRegionLine> generated;
    try (BufferedReader bufGen = new BufferedReader(generatedReader)) {
      generated = readCnvRegionLines(bufGen);
      if (mCollapse) {
        collapseCnvRegionLines(generated);
      }
    }

    final RegionPosition genPos = new RegionPosition(generated);
    final RegionPosition detPos = new RegionPosition(detected);
    while (detPos.valid() || genPos.valid()) {
      if (detPos.valid() && genPos.valid()) {
        if (closeEnough(genPos.getPos(), detPos.getPos())) {
          correctBreakpoint(genPos.getSequence(), genPos, detPos);
          genPos.next();
          detPos.next();
        } else if (genPos.getPos() < detPos.getPos()) {
          falseNegative(genPos.getSequence(), genPos);
          genPos.next();
        } else if (genPos.getPos() > detPos.getPos()) {
          falsePositive(genPos.getSequence(), detPos);
          detPos.next();
        }
      } else if (genPos.valid()) {
        falseNegative(genPos.getSequence(), genPos);
        genPos.next();
      } else {
        falsePositive(detPos.getSequence(), detPos);
        detPos.next();
      }
    }
    out.write("#RUN-ID\t" + CommandLine.getRunId() + StringUtils.LS);
    out.write("#true_positives\tfalse_positives\tfalse_negatives" + StringUtils.LS);
    out.write(mNumCorrect + TB + mNumFalsePositive + TB + mNumFalseNegative + StringUtils.LS);
  }

  private void correctBreakpoint(final String name, final RegionPosition genPos, final RegionPosition detPos) throws IOException {
    if (mCorrectOut != null) {
      mCorrectOut.write("correct: " + name + TB + genPos.getPos() + " ~ " + detPos.getPos() + "    means gen-det " + genPos.mRegion.mMean
          + " - " + detPos.mRegion.mMean + StringUtils.LS);
    }
    mNumCorrect++;
  }
  private void falseNegative(final String name, final RegionPosition genPos) throws IOException {
    if (mFalseNegativeOut != null) {
      mFalseNegativeOut.write("false negative: " + name + TB + genPos.getPos() + "    mean " + genPos.mRegion.mMean + StringUtils.LS);
    }
    mNumFalseNegative++;
  }
  private void falsePositive(final String name, final RegionPosition detPos) throws IOException {
    if (mFalsePositiveOut != null) {
      mFalsePositiveOut.write("false positive: " + name + TB + detPos.getPos() + "    mean " + detPos.mRegion.mMean + StringUtils.LS);
    }
    mNumFalsePositive++;
  }

  boolean closeEnough(final int genPos, final int detPos) {
    return Math.abs(genPos - detPos) <= mThreshold;
  }
}
