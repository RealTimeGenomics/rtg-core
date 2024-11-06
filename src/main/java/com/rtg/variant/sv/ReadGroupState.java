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
    if (!SamUtils.uniquelyMapped(rec)) {
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
