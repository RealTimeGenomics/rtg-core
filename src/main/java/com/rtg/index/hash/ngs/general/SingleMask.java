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
package com.rtg.index.hash.ngs.general;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds static information about all the bits in a single mask.
 * The information about each sub-mask is held in <code>Skel</code> objects.
 */
public final class SingleMask extends IntegralAbstract {
  /*
   * The number of chunks to be used &ge; number of skeletons to be stored
   * (two chunks without a gap between are stored as one skeleton).
   */
  private final int mChunks;

  /* number of bits to be set on summed over all chunks - same as window length. */
  private final int mWindowLength;

  private final Skel[] mSkeletons;

  private final int[] mGaps;

  /** Maximum allowed number of indels. */
  private final int mIndels;

  /** Maximum allowed length of an indel. */
  private final int mIndelLength;

  private boolean mFrozen = false;

  /** Number of sub-masks (&gt;= 1) */
  private int mSize = 0;

  /**
   * @return total number of bits in complete mask (summed over all sub-masks)
   */
  public int windowLength() {
    return mWindowLength;
  }

  /**
   * Get sub-skeleton.
   * @param i selects the sub-skeleton.
   * @return the sub-selection.
   */
  public Skel subSkeleton(final int i) {
    return mSkeletons[i];
  }

  /**
   * Get length of gap between <code>(i-1)'th and i'th</code> sub-mask.
   * For i == 0 assumes a fictitious sub-mask that ends at position -1.
   * @param i which sub-mask (&gt;= 0).
   * @return the length of the gap (&gt; 0 except for i==0 when &gt;= 0)
   */
  public int gaps(final int i) {
    return mGaps[i];
  }

  /**
   * @return the maximum number of allowed indels (&gt;= 0).
   */
  public int indels() {
    return mIndels;
  }

  /**
   * @return maximum allowed length of an indel (&gt;= 1).
   */
  public int indelLength() {
    return mIndelLength;
  }

  /**
   * @return the number of sub-masks (&gt;= 1 if frozen).
   */
  public int size() {
    return mSize;
  }

  /**
   * @param chunks number of chunks used (&ge; the number of skeletons).
   * @param windowLength number of bits to be set on summed over all all chunks - same as window length.
   * @param indels number of indels which are guaranteed to be found.
   * @param indelLength indel length
   */
  public SingleMask(final int chunks, final int windowLength, final int indels, final int indelLength) {
    mChunks = chunks;
    mWindowLength = windowLength;
    mSkeletons = new Skel[mChunks];
    mGaps = new int[mChunks];
    mIndels = indels;
    mIndelLength = indelLength;
    integrity();
  }

  /**
   * Add another skeleton.
   * @param sk skeleton to be added.
   */
  public void add(final Skel sk) {
    if (mFrozen) {
      throw new RuntimeException("adding to frozen MaskSel");
    }
    mSkeletons[mSize] = sk;
    final int gap;
    if (mSize == 0) {
      gap = sk.position() - sk.length() + 1;
      assert gap >= 0;
    } else {
      final Skel prev = mSkeletons[mSize - 1];
      gap = sk.position() - prev.position() - sk.length();
    }
    assert gap >= 0;
    mGaps[mSize] = gap;
    ++mSize;
  }

  /**
   * Prevent any more skeletons being added.
   */
  public void freeze() {
    if (mFrozen) {
      throw new RuntimeException();
    }
    mFrozen = true;
  }

  public boolean isFrozen() {
    return mFrozen;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Mask skeleton:").append(mFrozen ? "frozen" : "liquid").append(" indels=").append(mIndels).append(" indelLength=").append(mIndelLength).append(LS);
    for (int i = 0; i < mSize; ++i) {
      final Skel sk = mSkeletons[i];
      sb.append(mGaps[i]).append("[").append(i).append("]").append(sk).append(LS);
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mGaps != null && mSkeletons != null && mGaps.length == mSkeletons.length);
    Exam.assertTrue(mChunks >= 0);
    Exam.assertTrue(mWindowLength >= 0 && mWindowLength <= 32);
    Exam.assertTrue(mSize >= 0 && mSize <= mSkeletons.length);

    int total = 0;
    int posn = -1;
    for (int i = 0; i < mSize; ++i) {
      final Skel sk = mSkeletons[i];
      total += sk.length();
      posn += sk.length() + mGaps[i];
      Exam.assertEquals(posn, sk.position());
      Exam.assertEquals(total, sk.finalPosition() + 1);
    }
    Exam.assertTrue(total <= mWindowLength);

    if (mFrozen) {
      Exam.assertTrue(mChunks >= mSize);
      Exam.assertTrue(mSize > 0);
      Exam.assertEquals(total, mWindowLength);
    }
    for (int i = 0; i < mSize - 1; ++i) {
      final int gap = mGaps[i];
      Exam.assertTrue(i == 0 ? gap >= 0 : gap > 0);
    }
    return true;
  }
}
