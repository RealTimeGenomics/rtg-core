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

package com.rtg.visualization;

import java.io.File;

import com.rtg.util.intervals.RegionRestriction;

/**
 */
class AviewParamsBuilder {

  private static final File[] EMPTY = new File[0];

  protected File[] mAlignments = EMPTY;
  protected File[] mUnmapped = EMPTY;
  protected File[] mTrackFiles = EMPTY;
  protected File mBaselineFile = null;
  protected File mReference = null;
  protected File[] mReads = EMPTY;

  protected boolean mDisplayDots = true;
  protected boolean mUseTerminalColor = true;
  protected boolean mColorBases = false;
  protected boolean mUseHtml = false;
  protected boolean mPrintCigars = false;
  protected boolean mPrintReadName = false;
  protected boolean mPrintMapQ = false;
  protected boolean mSortReads = false;
  protected boolean mSortReadGroup = false;
  protected boolean mSortSample = false;
  protected boolean mPrintReadGroup = false;
  protected boolean mPrintSample = false;
  protected boolean mPrintMatePosition = false;
  protected boolean mUnflattenCgi = false;
  protected boolean mShowSoftClippedBases = false;
  protected String[] mSamples = null;

  protected RegionRestriction mRegion = new RegionRestriction("", RegionRestriction.MISSING, RegionRestriction.MISSING);
  protected int mRegionPadding = -1;
  protected int mHeaderLineRepeat = 5;
  protected int mMaxMatedAlignmentScore = Integer.MAX_VALUE;
  protected int mMaxUnmatedAlignmentScore = Integer.MAX_VALUE;
  protected int mMaxIhScore = Integer.MAX_VALUE;
  protected int mMinMapQ = 0;
  protected int mMappingTolerance = 0;
  protected int mProjectTrackId = 0;


  AviewParams create() {
    return new AviewParams(this);
  }

  AviewParamsBuilder alignments(File[] alignments) {
    mAlignments = alignments;
    return this;
  }

  AviewParamsBuilder samples(String[] samples) {
    mSamples = samples;
    return this;
  }

  AviewParamsBuilder unmapped(File[] unmapped) {
    mUnmapped = unmapped;
    return this;
  }

  AviewParamsBuilder reads(File[] reads) {
    mReads = reads;
    return this;
  }

  AviewParamsBuilder baselineFile(File baselineFile) {
    mBaselineFile = baselineFile;
    return this;
  }

  AviewParamsBuilder trackFiles(File[] trackFiles) {
    mTrackFiles = trackFiles;
    return this;
  }

  AviewParamsBuilder reference(File reference) {
    mReference = reference;
    return this;
  }

  AviewParamsBuilder displayDots(boolean displayDots) {
    mDisplayDots = displayDots;
    return this;
  }

  AviewParamsBuilder useTerminalColor(boolean colorTerminal) {
    mUseTerminalColor = colorTerminal;
    return this;
  }

  AviewParamsBuilder colorBases(boolean colorBases) {
    mColorBases = colorBases;
    return this;
  }

  AviewParamsBuilder unflattenCgi(boolean unflatten) {
    mUnflattenCgi = unflatten;
    return this;
  }

  AviewParamsBuilder projectTrackId(int trackId) {
    mProjectTrackId = trackId;
    return this;
  }

  AviewParamsBuilder useHtml(boolean html) {
    mUseHtml = html;
    return this;
  }

  AviewParamsBuilder printCigars(boolean printCigars) {
    mPrintCigars = printCigars;
    return this;
  }

  AviewParamsBuilder printReadName(boolean printReadName) {
    mPrintReadName = printReadName;
    return this;
  }

  AviewParamsBuilder region(String region) {
    mRegion = new RegionRestriction(region);
    return this;
  }

  AviewParamsBuilder regionPadding(int padding) {
    mRegionPadding = padding;
    return this;
  }

  AviewParamsBuilder headerLineRepeat(int headerLineRepeat) {
    mHeaderLineRepeat = headerLineRepeat;
    return this;
  }

  AviewParamsBuilder maxMatedAlignmentScore(int maxAlignmentScore) {
    mMaxMatedAlignmentScore = maxAlignmentScore;
    return this;
  }

  AviewParamsBuilder maxUnmatedAlignmentScore(int maxAlignmentScore) {
    mMaxUnmatedAlignmentScore = maxAlignmentScore;
    return this;
  }

  AviewParamsBuilder maxIhScore(int maxIhScore) {
    mMaxIhScore = maxIhScore;
    return this;
  }

  AviewParamsBuilder minMapQ(int minMapq) {
    mMinMapQ = minMapq;
    return this;
  }

  AviewParamsBuilder sortReads(boolean sortReads) {
    mSortReads = sortReads;
    return this;
  }

  AviewParamsBuilder sortReadGroup(boolean sortReadGroup) {
    mSortReadGroup = sortReadGroup;
    return this;
  }

  AviewParamsBuilder sortSample(boolean sortSample) {
    mSortSample = sortSample;
    return this;
  }

  AviewParamsBuilder printReadGroup(boolean printReadGroup) {
    mPrintReadGroup = printReadGroup;
    return this;
  }

  AviewParamsBuilder printSample(boolean print) {
    mPrintSample = print;
    return this;
  }

  AviewParamsBuilder printMapQ(boolean printMapQ) {
    mPrintMapQ = printMapQ;
    return this;
  }

  AviewParamsBuilder printMatePosition(boolean printMatePosition) {
    mPrintMatePosition = printMatePosition;
    return this;
  }

  AviewParamsBuilder mappingTolerance(int tolerance) {
    mMappingTolerance = tolerance;
    return this;
  }

  AviewParamsBuilder showSoftClippedBases(boolean val) {
    mShowSoftClippedBases = val;
    return this;
  }
}
