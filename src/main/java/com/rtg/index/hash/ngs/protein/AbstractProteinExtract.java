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
package com.rtg.index.hash.ngs.protein;

import java.io.IOException;

import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.index.hash.ngs.general.Skel;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Deals with extracting masks.
 * Indel tweaking is done dynamically in the <code>maskIndel()</code> method.
 */
public abstract class AbstractProteinExtract extends IntegralAbstract {

  protected final SingleMask mSingleMask;

  AbstractProteinExtract(final SingleMask singleMask) {
    mSingleMask = singleMask;
  }

  /**
   * Called for each mask extraction.
   * Needed to accomodate the multiple masks that can result from allowance for indels.
   * @param bits final extracted bits using mask.
   * @throws IOException If an I/O error occurs
   */
  protected abstract void masked(final long bits) throws IOException;

  /**
   * Extract the bits from a single word.
   * @param v word bits are to be extracted from.
   * @return the extracted bits (right shifted aligned to position 0).
   */
  long mask(final long v) {
    assert mSingleMask.isFrozen();
    long res = 0;
    for (int i = 0; i < mSingleMask.size(); ++i) {
      final Skel sk = mSingleMask.subSkeleton(i);
      res = res | sk.mask(v);
    }
    return res;
  }

  /**
   * Extract the bits from an array of words and pack them into a single word.
   * @param vs the words bits are extracted from (most significant bit first).
   * @throws IOException If an I/O error occurs
   */
  void mask(long[] vs) throws IOException {
    assert mSingleMask.isFrozen();
    assert vs.length * mSingleMask.windowLength() <= 64;
    //System.err.println(mSingleMask.mWindowLength);
    long res = 0L;
    for (long v : vs) {
      res = res << mSingleMask.windowLength() | mask(v);
    }
    masked(res);
  }

  /**
   * Extract the bits from an array of words and pack the result into a single word.
   * Includes indel tweaking.
   *
   * @param vs the words bits are extracted from (most significant bit first).
   * @throws IOException If an I/O error occurs
   */
  //TODO this code could probably do with some micro-optimization
  void maskIndel(final long[] vs) throws IOException {
    assert mSingleMask.isFrozen();
    ma(vs, 0, 0, 0, 0, 0, 0, mSingleMask.indels(), 0);
  }

  /*
   * @param rs rs[i] is the accumulated bits from <code>vs[i]</code> so far
   * @param i index of sub-mask
   * @param indels number of indels remaining to be used
   * @param tweak tweak parameter
   * @throws IOException If an I/O error occurs
   */
  private void ma(final long[] vs, final long a, final long b, final long c, final long d, final long e, int i, int indels, int tweak) throws IOException {
    assert indels >= 0;
    assert tweak <= mSingleMask.indels() * mSingleMask.indelLength() && -mSingleMask.indels()  * mSingleMask.indelLength() <= tweak;
    assert 0 <= i && i < mSingleMask.size();

    final Skel sk = mSingleMask.subSkeleton(i);
    final long aa = a | sk.mask(vs[0], tweak);
    final long bb = b | sk.mask(vs[1], tweak);
    final long cc = c | sk.mask(vs[2], tweak);
    final long dd = d | sk.mask(vs[3], tweak);
    final long ee = e | sk.mask(vs[4], tweak);

    final int j = i + 1;
    if (j == mSingleMask.size()) {
      //System.err.println("done i=" + i);
      //combine the sets of extracted bits into a single word
      long res = 0L;
      res = (res << mSingleMask.windowLength()) | aa;
      res = (res << mSingleMask.windowLength()) | bb;
      res = (res << mSingleMask.windowLength()) | cc;
      res = (res << mSingleMask.windowLength()) | dd;
      res = (res << mSingleMask.windowLength()) | ee;
      masked(res);
      return;
    }
    final int indelLength = mSingleMask.indelLength();
    final int indGap = indels * indelLength;
    final int min = Math.min(indGap, mSingleMask.gaps(i + 1));
    //System.err.println("indels=" + indels + " indelLength=" + indelLength + " indGap=" + indGap + " min=" + min);
    for (int t = -min; t < 0; ++t) {
      final int s = (-t + indelLength - 1) / indelLength;
      //System.err.println("neg i=" + i + " t=" + t + " s=" + s);
      ma(vs, aa, bb, cc, dd, ee, j, indels - s, tweak + t);
    }
    for (int t = 0; t <= indGap; ++t) {
      final int s = (t + indelLength - 1) / indelLength;
      //System.err.println("pos i=" + i + " t=" + t + " s=" + s);
      ma(vs, aa, bb, cc, dd, ee, j, indels - s, tweak + t);
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ProteinExtract");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSingleMask != null);
    return true;
  }
}
