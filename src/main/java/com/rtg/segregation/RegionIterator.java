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

package com.rtg.segregation;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.util.Pair;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 *
 */
@TestClass("com.rtg.segregation.SegregationVcfSearchTest")
class RegionIterator extends IntegralAbstract {

  private final PrintStream mOut;
  private final BedWriter mBed;
  private final String mSeq;
  private final int mLength;
  private int mStart = -1;
  private int mEnd = -1;
  private SearchType mStartType = null;
  private PatternArray mPa = null;
  private int mFaLevel = 0;
  private int mMoLevel = 0;
  private final Map<String, Integer> mChrLengths;
  private final Map<Pair<Sex, String>, ReferenceSequence> mPloidyMap;

  RegionIterator(PrintStream out, BedWriter bed, String seq, int length, Map<String, Integer> chrLengths, Map<Pair<Sex, String>, ReferenceSequence> ploidyMap) {
    mOut = out;
    mBed = bed;
    mSeq = seq;
    mLength = length;
    mChrLengths = chrLengths;
    mPloidyMap = ploidyMap;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mOut);
    Exam.assertNotNull(mBed);
    Exam.assertNotNull(mSeq);
    Exam.assertTrue(mLength >= 1);
    Exam.assertTrue(0 <= mFaLevel && mFaLevel <= mLength);
    Exam.assertTrue(0 <= mMoLevel && mMoLevel <= mLength);
    if (mStartType == null) {
      Exam.assertEquals(-1, mStart);
      Exam.assertEquals(-1, mEnd);
      Exam.assertTrue(mPa == null);
    } else {
      Exam.assertTrue(0 <= mStart && mStart < mEnd);
      Exam.assertNotNull(mPa);
    }
    return true;
  }

  int minDiff(final int a, final int b) {
    final int c = mLength - b;
    final int min0 = Math.abs(a - b);
    final int min1 = Math.abs(a - c);
    if (min0 <= min1) {
      return b;
    } else {
      return c;
    }
  }

  void iterate(final Iterator<SearchContainer> it) throws IOException {
    //System.err.println(toString());
    while (it.hasNext()) {
      assert integrity();
      final SearchContainer sc = it.next();
      final SearchType type = sc.type();
      final SegregationBlock block = sc.block();
      //System.err.println(type.toString());
      //System.err.println(block.toString());
      mOut.println(type.code() + " " + block + " " + sc.pattern().faString() + " " + sc.pattern().moString());
      switch (type) {
        case OK:
          mEnd = block.end();
          mPa = sc.pattern();
          mFaLevel = minDiff(mFaLevel, mPa.faCount());
          mMoLevel = minDiff(mMoLevel, mPa.moCount());
          if (mStartType == null) {
            //An OK instead of a New at the start (it was just too hard to get it right).
            mStartType = SearchType.New;
            //mStart = block.start() - 1;
            mStart = 0;
            mEnd = block.end();
          }
          break;
        case Error:
          //do nothing
          break;
        case XO:
          final CrossOver xo = sc.xo();
          //TODO this output is redundant and can be got rid of - but is currently used by some of the downstream scripts
          mOut.println("XO " + mSeq + " " + block.start() + " " + xo);
          out(block.start() - 1);
          mStart = mEnd;
          mStartType = SearchType.XO;
          mEnd = block.end();
          mPa = sc.pattern();
          final int[] sl = xo.searchLevel(mFaLevel, mMoLevel);
          mFaLevel = sl[0];
          mMoLevel = sl[1];
          break;
        case New:
          out(mEnd);
          mStart = mStartType == null ? 0 : block.start() - 1;
          mStartType = SearchType.New;
          mEnd = block.end();
          mPa = sc.pattern();
          mFaLevel = minDiff(mFaLevel, mPa.faCount());
          mMoLevel = minDiff(mMoLevel, mPa.moCount());
          break;
        default:
          throw new RuntimeException();
      }
      //System.err.println(toString());
      assert integrity();
    }
    out(mChrLengths.get(mSeq));
  }

  void out(int end) throws IOException {
    if (mPa == null || !mPa.isUniqueAndValid(isXLike())) {
      return;
    }
    //bed is zero based and vcf is 1 based
    mBed.write(new BedRecord(mSeq, mStart, end, mPa.faString(), mPa.moString(), mStartType.bedCode(), String.valueOf(mFaLevel), String.valueOf(mMoLevel)));
  }

  private boolean isXLike() {
    final Ploidy female = mPloidyMap.get(new Pair<>(Sex.FEMALE, mSeq)).effectivePloidy(mStart);
    final Ploidy male = mPloidyMap.get(new Pair<>(Sex.MALE, mSeq)).effectivePloidy(mStart);
    return male == Ploidy.HAPLOID && female == Ploidy.DIPLOID;
  }
}
