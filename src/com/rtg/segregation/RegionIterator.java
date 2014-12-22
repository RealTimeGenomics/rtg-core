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
      mOut.println(type.code() + " " + block.toString() + " " + sc.pattern().faString() + " " + sc.pattern().moString());
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
          mOut.println("XO " + mSeq + " " + block.start() + " " + xo.toString());
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
    mBed.write(new BedRecord(mSeq, mStart, end, mPa.faString(), mPa.moString(), mStartType.bedCode(), "" + mFaLevel, "" + mMoLevel));
  }

  private boolean isXLike() {
    final Ploidy female = mPloidyMap.get(new Pair<>(Sex.FEMALE, mSeq)).effectivePloidy(mStart);
    final Ploidy male = mPloidyMap.get(new Pair<>(Sex.MALE, mSeq)).effectivePloidy(mStart);
    return male == Ploidy.HAPLOID && female == Ploidy.DIPLOID;
  }
}
