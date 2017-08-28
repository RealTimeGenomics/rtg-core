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
package com.rtg.variant.cnv.region;

import java.util.ArrayList;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.StringUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
final class ComplexCnvRegion extends AbstractCnvRegion implements Integrity {

  private final AbstractCnvRegion mFirst;

  private final AbstractCnvRegion[] mMiddle;

  private final int mLengthMid;

  private final int mSpanMid;

  private final AbstractCnvRegion mLast;

  ComplexCnvRegion(final int start, final int end, final SortedSet<AbstractCnvRegion> regions) {
    super(start, end);
    final int size = regions.size();
    if (size <= 1) {
      throw new IllegalArgumentException();
    }
    final AbstractCnvRegion[] rega = regions.toArray(new AbstractCnvRegion[size]);
    mLengthMid = 2 * size;
    mSpanMid = getEnd() - getStart() + 1;
    mFirst = rega[0];
    mLast = rega[size - 1];
    if (size == 2) {
      mMiddle = null;
    } else {
      final ArrayList<SortedSet<AbstractCnvRegion>> sets = new ArrayList<>();
      for (int i = 0; i < mLengthMid; ++i) {
        sets.add(new TreeSet<AbstractCnvRegion>());
      }
      for (int i = 1; i < size - 1; ++i) {
        final AbstractCnvRegion reg = rega[i];
        final int bs = box(Math.max(reg.getStart(), getStart()));
        final int be = box(Math.min(reg.getEnd(), getEnd()));
        for (int j = bs; j <= be; ++j) {
          sets.get(j).add(reg);
        }
      }
      mMiddle = new AbstractCnvRegion[mLengthMid];
      for (int i = 0; i < mLengthMid; ++i) {

        final SortedSet<AbstractCnvRegion> set = sets.get(i);
        final int ss = set.size();
        if (ss == 0) {
          continue;
        }
        if (ss == 1) {
          mMiddle[i] = set.first();
          continue;
        }
        mMiddle[i] = new ComplexCnvRegion(index(i), index(i + 1) - 1, set);
      }
    }
  }

  int box(final int index) {
    assert index >= getStart() && index <= getEnd();
    return (int) ((index - getStart()) * (long) mLengthMid / mSpanMid);
  }

  int index(final int box) {
    assert box >= 0 && box <= mLengthMid;
    return getStart() + (int) (((long) box * mSpanMid + mLengthMid - 1) / mLengthMid);
  }

  @Override
  public boolean isInRegion(final int index) {
    if (index < getStart() || index > getEnd()) {
      return false;
    }
    if (mFirst.isInRegion(index)) {
      return true;
    }
    if (mLast.isInRegion(index)) {
      return true;
    }
    if (mMiddle == null) {
      return false;
    }
    final int box = box(index);
    if (mMiddle[box] == null) {
      return false;
    }
    return mMiddle[box].isInRegion(index);
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("ComplexRegion start=").append(getStart()).append(" end=").append(getEnd()).append(" length=").append(mLengthMid).append(StringUtils.LS);
    sb.append("First ");
    sb.append(mFirst);
    sb.append(StringUtils.LS);
    if (mMiddle != null) {
      for (int i = 0; i < mMiddle.length; ++i) {
        if (mMiddle[i] != null) {
          sb.append("[").append(i).append("]").append(StringUtils.LS);
          sb.append(mMiddle[i]);
          sb.append(StringUtils.LS);
        }
      }
    }
    sb.append("Last ");
    sb.append(mLast);
    sb.append(StringUtils.LS);
    return sb.toString();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    Exam.globalIntegrity(mFirst);
    Exam.globalIntegrity(mLast);
    if (mMiddle != null) {
      for (final AbstractCnvRegion element : mMiddle) {
        if (element != null) {
          Exam.globalIntegrity(element);
          checkSubRegion(element);
        }
      }
    }
    return true;
  }

  void checkSubRegion(final AbstractCnvRegion r) {
    Exam.assertNotNull(r);
    Exam.assertFalse(r.getEnd() < getStart());
    Exam.assertFalse(r.getStart() > getEnd());
  }

  @Override
  public boolean integrity() {
    checkSubRegion(mFirst);
    checkSubRegion(mLast);
    Exam.assertEquals(mSpanMid, getEnd() - getStart() + 1);
    Exam.assertTrue(mLengthMid >= 4);
    if (mMiddle != null) {
      Exam.assertEquals(mLengthMid, mMiddle.length);
    }
    return true;
  }

}
