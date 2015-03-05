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

package com.rtg.variant.sv;


import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMRecord;

/**
 * Encapsulates per read group state maintained during structural variant calling.
 *
 */
public class ReadGroupState {

  private final ReadGroupStats mStats;
  private final MachineType mMachine;
  private final MachineOrientation mOrientation;

  private final SamCounts mNotPaired;

  private final SamCounts mUnmatedLeftArm;
  private final SamCounts mProperLeftArm;
  private final SamCounts mDiscordantLeftArm;

  private final SamCounts mUnmatedRightArm;
  private final SamCounts mProperRightArm;
  private final SamCounts mDiscordantRightArm;
  private final SamCounts mUniqueCoverage;
  private final SamCounts mAmbiguousCoverage;

  private final int mLo;
  private final int mHi;


  ReadGroupState(ReadGroupStats stats, MachineType machineType, Corrections corrections) {
    mStats = stats;
    mMachine = machineType;
    mOrientation = mMachine.orientation();
    mLo = stats.lo();
    mHi = stats.hi();
    mNotPaired = new CumulativeSamCounts(0, corrections);
    mUnmatedLeftArm = new CumulativeSamCounts(0, corrections);
    mProperLeftArm = new CumulativeSamCounts(0, corrections);
    mDiscordantLeftArm = new CumulativeSamCounts(0, corrections);
    mUnmatedRightArm = new CumulativeSamCounts(0, corrections);
    mProperRightArm = new CumulativeSamCounts(0, corrections);
    mDiscordantRightArm = new CumulativeSamCounts(0, corrections);
    mUniqueCoverage = new CumulativeSamCounts(0, corrections);
    mAmbiguousCoverage = new CumulativeSamCounts(0, corrections);
  }

  //This version of constructor only used in testing
  ReadGroupState(ReadGroupStats stats, AllCounts counts) {
    mStats = stats;
    mMachine = MachineType.ILLUMINA_PE;
    mOrientation = mMachine.orientation();
    mLo = stats.lo();
    mHi = stats.hi();
    Diagnostic.developerLog("Window size is " + mLo + ":" + mHi);
    mNotPaired = counts.unpaired();
    mUnmatedLeftArm = counts.unmatedLeft();
    mProperLeftArm = counts.properLeft();
    mDiscordantLeftArm = counts.discordantLeft();
    mUnmatedRightArm = counts.unmatedRight();
    mProperRightArm = counts.properRight();
    mDiscordantRightArm = counts.discordantRight();
    mUniqueCoverage = new CumulativeSamCounts(0, null);
    mAmbiguousCoverage = new CumulativeSamCounts(0, null);
  }

  ReadGroupStats stats() {
    return mStats;
  }

  int lo() {
    return mLo;
  }

  int hi() {
    return mHi;
  }

  MachineType machine() {
    return mMachine;
  }

  MachineOrientation orientation() {
    return mOrientation;
  }

  SamCounts notPaired() {
    return mNotPaired;
  }

  SamCounts unique() {
    return mUniqueCoverage;
  }

  SamCounts ambiguous() {
    return mAmbiguousCoverage;
  }

  SamCounts properLeftArm(boolean reverse) {
    return reverse ? mProperRightArm : mProperLeftArm;
  }

  SamCounts discordantLeftArm(boolean reverse) {
    return reverse ? mDiscordantRightArm : mDiscordantLeftArm;
  }

  SamCounts unmatedLeftArm(boolean reverse) {
    return reverse ? mUnmatedRightArm : mUnmatedLeftArm;
  }

  SamCounts properRightArm(boolean reverse) {
    return reverse ? mProperLeftArm : mProperRightArm;
  }

  SamCounts discordantRightArm(boolean reverse) {
    return reverse ? mDiscordantLeftArm : mDiscordantRightArm;
  }

  SamCounts unmatedRightArm(boolean reverse) {
    return reverse ? mUnmatedLeftArm : mUnmatedRightArm;
  }

  void reset(int length, int bufferSize) {
    mNotPaired.reset(length, bufferSize);
    mUnmatedLeftArm.reset(length, bufferSize);
    mProperLeftArm.reset(length, bufferSize);
    mDiscordantLeftArm.reset(length, bufferSize);
    mUnmatedRightArm.reset(length, bufferSize);
    mProperRightArm.reset(length, bufferSize);
    mDiscordantRightArm.reset(length, bufferSize);
    mUniqueCoverage.reset(length, bufferSize);
    mAmbiguousCoverage.reset(length, bufferSize);
  }

  void flushTo(int offset) {
    mNotPaired.flushTo(offset);
    mUnmatedLeftArm.flushTo(offset);
    mProperLeftArm.flushTo(offset);
    mDiscordantLeftArm.flushTo(offset);
    mUnmatedRightArm.flushTo(offset);
    mProperRightArm.flushTo(offset);
    mDiscordantRightArm.flushTo(offset);
    mUniqueCoverage.flushTo(offset);
    mAmbiguousCoverage.flushTo(offset);
  }

  boolean update(SAMRecord rec) {
    final int start = rec.getAlignmentStart() - 1; // zero-based, inclusive
    Integer nh = SamUtils.getNHOrIH(rec);
    if (nh == null) {
      //already should have warned in svprep
      nh = 1;
    }
    if (nh > 1) {
      mAmbiguousCoverage.increment(start);
      return true;
    }
    mUniqueCoverage.increment(start);
    final boolean readPaired = rec.getReadPairedFlag();
    if (!readPaired) {
      mNotPaired.increment(start);
      return true;
    }
    final boolean mateUnmapped = rec.getMateUnmappedFlag();
    final boolean first = mOrientation.firstOnTemplate(rec);
    if (mateUnmapped) {
      (first ? mUnmatedLeftArm : mUnmatedRightArm).increment(start);
      return true;
    }
    if (!mOrientation.hasValidMate(rec)) {
      (first ? mDiscordantLeftArm : mDiscordantRightArm).increment(start);
      return true;
    }
    final int insertSize = Math.abs(rec.getInferredInsertSize());
    if (inNormalDistribution(insertSize, mStats)) {
      (first ? mProperLeftArm : mProperRightArm).increment(start);
      return true;
    }
    (first ? mDiscordantLeftArm : mDiscordantRightArm).increment(start);
    return true;
  }

  private static boolean inNormalDistribution(int insertSize, ReadGroupStats stats) {
    final double x = Math.abs((insertSize - stats.fragmentMean()) / stats.fragmentStdDev());
    if (x < 3.0) {
      return true;
    }
    return false;
  }
}
