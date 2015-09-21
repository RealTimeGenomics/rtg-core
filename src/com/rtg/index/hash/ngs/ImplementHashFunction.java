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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.index.hash.PrimeUtils;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Tailored for use with <code>Ngs</code> processing.
 * Allows split windows to be extracted.
 */
public abstract class ImplementHashFunction extends IntegralAbstract implements NgsHashFunction, Cloneable {

  /** Number of bits in an int */
  private static final int INT_BITS = 32;

  /** Number of bits in a long */
  protected static final int LONG_BITS = 64;

  protected final int mReadLength;

  /** Number of bits to right shift when using reverse complement buffers for reads. */
  protected final int mReadLengthReverse;

  /** Number of bits to right shift when using reverse complement buffers for template (differs from <code>mReadLengthReverse</code> for CG. */
  protected final int mTemplateLengthReverse;

  /** Number of codes to be included in a hash. */
  private final int mWindowSize;

  private final int mWindowBits;

  private final long mPrime;

  protected final MemScore mMemScore;

  protected long[] mReadSequencesF1 = null;
  protected long[] mReadSequencesF2 = null;

  /** Mask of the recent bits covering exactly the window. */
  private final long mMask;

  protected ReadCall mReadCall;

  protected TemplateCall mTemplateCall;

  //below here is dynamic

  /**
   * High order bits of window.
   */
  protected long mValuesF0;

  /**
   * Low order bits of window.
   */
  protected long mValuesF1;

  /**
   * High order bits of window.
   */
  protected long mValuesR0;

  /**
   * Low order bits of window.
   */
  protected long mValuesR1;

  protected int mSoFar = 0;

  /**
   * @param readLength number of nucleotides in a complete read.
   * @param windowSize number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public ImplementHashFunction(final int readLength, final int windowSize, final ReadCall readCall, final TemplateCall templateCall) {
    this(readLength, windowSize, readCall, templateCall, LONG_BITS - readLength);
  }

  /**
   * @param readLength number of nucleotides in a complete read.
   * @param windowSize number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   * @param readLengthReverse how much to shift right when in reverse complement. Needs to be specified explicitly for CG masks.
   */
  public ImplementHashFunction(final int readLength, final int windowSize, final ReadCall readCall, final TemplateCall templateCall, final int readLengthReverse) {
    //System.err.println("ImplementHashFunction readLength=" + readLength + " windowSize=" + windowSize);
    mReadLength = readLength;
    mReadLengthReverse = LONG_BITS - readLength;
    mTemplateLengthReverse = readLengthReverse;
    mWindowSize = windowSize;
    mReadCall = readCall;
    mTemplateCall = templateCall;
    mWindowBits = mWindowSize * 2;
    mMask = mWindowBits == 64 ? -1L : (1L << mWindowBits) - 1;
    mPrime = PrimeUtils.prime(mWindowBits);
    mMemScore = new MemScore(readLength());
    resetPrivate();
  }

  @Override
  public NgsHashFunction threadClone(final HashingRegion region) throws IOException {
    try {
      final ImplementHashFunction clone = (ImplementHashFunction) clone();
      clone.mTemplateCall = this.mTemplateCall.threadClone(region);
      clone.setHashFunction();
      return clone;
    } catch (final CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public void threadFinish() throws IOException {
    mTemplateCall.threadFinish();
    mReadSequencesF1 = null;
    mReadSequencesF2 = null;
  }

  @Override
  public void templateBidirectional(final int endPosition) throws IOException {
    templateForward(endPosition);
    templateReverse(endPosition);
  }

  @Override
  public void templateForward(final int endPosition) throws IOException {
    if (mSoFar >= mReadLength) { // XXXLen This is an assumption that the template is spanned by mReadLength bases, not true for inserts or for CG reads.
      mTemplateCall.setReverse(false);
      templateAll(endPosition, mValuesF0, mValuesF1);
    }
  }

  @Override
  public void templateReverse(final int endPosition) throws IOException {
    // Only start when reverse is full, adjust end position
    if (mSoFar >= 64) {
      mTemplateCall.setReverse(true);
      templateAll(endPosition - mTemplateLengthReverse, mValuesR0, mValuesR1);
    }
  }

  @Override
  public void endSequence() throws IOException {
    mTemplateCall.endSequence();
  }

  /**
   * Process all windows for the specified position in the template.
   * @param endPosition last position in the current template window.
   * @param v0 low order bits for each nucleotide in the current window
   * @param v1 high order bits for each nucleotide in the current window
   * @throws IOException If an I/O error occurs
   */
  public abstract void templateAll(int endPosition, long v0, long v1) throws IOException;


  @Override
  public void readAll(final int readId, final boolean reverse) throws IOException {
    //System.err.println("readAll readId=" + readId + " reverse=" + reverse);
    if (reverse) {
      readAll(readId, mValuesR0 >>> mReadLengthReverse, mValuesR1 >>> mReadLengthReverse);
    } else {
      readAll(readId, mValuesF0, mValuesF1);
    }
  }


  /**
   * Process all windows for the specified read.
   * @param readId number of the read (&gt;=0).
   * @param v0 low order bits for each nucleotide in the current window
   * @param v1 high order bits for each nucleotide in the current window
   * @throws IOException If an I/O error occurs
   */
  public abstract void readAll(int readId, long v0, long v1) throws IOException;

  /**
   * Tell the template about the hash function.
   * Done here to ensure it is after the subclass
   * constructors are complete.
   */
  protected void setHashFunction() {
    mTemplateCall.setHashFunction(this);
  }

  @Override
  public void hashStep(final byte code) {
    assert code >= 0 && code <= 3;
    final int c1 = code >> 1;
    final int c2 = code & 1;
    mValuesF0 = (mValuesF0 << 1) | c1;
    mValuesF1 = (mValuesF1 << 1) | c2;

    final int cr = 3 - code;
    final long cr1 = cr >> 1;
    final long cr2 = cr & 1;
    mValuesR0 = mValuesR0 >>> 1;
    mValuesR0 |= cr1 << 63;
    mValuesR1 = mValuesR1 >>> 1;
    mValuesR1 |= cr2 << 63;

    mSoFar++;
    //System.err.println("hashStep code=" + code);
    //System.err.println("F0=" + com.rtg.util.Utils.toBitsSep(mValuesF0) + "\nF1=" + com.rtg.util.Utils.toBitsSep(mValuesF1));
    //System.err.println("R0=" + com.rtg.util.Utils.toBitsSep(mValuesR0) + "\nR1=" + com.rtg.util.Utils.toBitsSep(mValuesR1));
  }

  @Override
  public void hashStep() {
    hashStep((byte) 0);
  }

  @Override
  public int readLength() {
    return mReadLength;
  }

  @Override
  public int windowSize() {
    return mWindowSize;
  }

  /** This is called from constructor, so must be private. */
  private void resetPrivate() {
    mSoFar = 0;
    mValuesF0 = 0;
    mValuesF1 = 0;
    mValuesR0 = -1L;
    mValuesR1 = -1L;
  }

  @Override
  public void reset() {
    resetPrivate();
  }

  /**
   * Hash the specified value.
   * Over-ridden in testing so can see underlying value.
   * @param x value to be hashed.
   * @return hashed value of x.
   */
  protected long hash(final long x) {
    return (x * mPrime) & mMask;
  }

  @Override
  public void setReadSequences(final long numberReads) {
    if (numberReads > Integer.MAX_VALUE || numberReads < 0) {
      throw new RuntimeException("numberReads=" + numberReads);
    }
    mReadSequencesF1 = new long[(int) numberReads];
    mReadSequencesF2 = new long[(int) numberReads];
  }

  @Override
  public void templateSet(final long name, final int length) {
    mTemplateCall.set(name, length);
  }

  @Override
  public void setValues(final int id2, final boolean reverse) {
    //System.err.println("setValues id2=" + id2 + " reverse=" + reverse);
    if (reverse) {
      mReadSequencesF1[id2] = mValuesR0 >>> mReadLengthReverse;
      mReadSequencesF2[id2] = mValuesR1 >>> mReadLengthReverse;
    } else {
      mReadSequencesF1[id2] = mValuesF0;
      mReadSequencesF2[id2] = mValuesF1;
    }
    //System.err.println("0=" + com.rtg.util.Utils.toBitsSep(mReadSequencesF[id2]));
    //System.err.println("1=" + com.rtg.util.Utils.toBitsSep(mReadSequencesF[id2 + 1]));
  }

  @Override
  public int fastScore(final int readId) {
    //System.err.println("fastScore readId=" + readId + " rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    if (mTemplateCall.isReverse()) {
      return mMemScore.fastScore(bit0, bit1, mValuesR0, mValuesR1);
    } else {
      return mMemScore.fastScore(bit0, bit1, mValuesF0, mValuesF1);
    }
  }

  @Override
  public int indelScore(final int readId) {
    //System.err.println("indelScore rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    if (mTemplateCall.isReverse()) {
      return mMemScore.indelScore(bit0, bit1, mValuesR0, mValuesR1);
    } else {
      return mMemScore.indelScore(bit0, bit1, mValuesF0, mValuesF1);
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ImplementHashFunction read length=").append(mReadLength).append(" window size=").append(mWindowSize);
  }

  /**
   * Check internal state.
   * Intended for unit-tests.
   * @param a sequence of high order bits of codes.
   * @param b sequence of low order bits of codes.
   * @param soFar number of valid <code>hashStep</code> calls since last reset.
   */
  void integrity(final long a, final long b, final int soFar) {
    //System.err.println(" a=" + com.rtg.util.Utils.toBitsSep(a) + "  b=" + com.rtg.util.Utils.toBitsSep(b));
    //System.err.println("f0=" + com.rtg.util.Utils.toBitsSep(mValuesF0) + " f1=" + com.rtg.util.Utils.toBitsSep(mValuesF1));
    //System.err.println("r0=" + com.rtg.util.Utils.toBitsSep(mValuesR0) + " r1=" + com.rtg.util.Utils.toBitsSep(mValuesR1));
    Exam.assertTrue(a == mValuesF0);
    Exam.assertTrue(b == mValuesF1);
    Exam.assertTrue(bitFlip(a) == mValuesR0);
    Exam.assertTrue(bitFlip(b) == mValuesR1);

    Exam.assertTrue("soFar=" + soFar + " mSoFar=" + mSoFar, soFar == mSoFar);
  }

  static long bitFlip(final long a) {
    long ra = a;
    long na = 0;
    for (int i = 0; i < 64; i++) {
      na = (na << 1) | (1 - (ra & 1));
      ra = ra >> 1;
    }
    return na;
  }

  @Override
  public void logStatistics() {
    mTemplateCall.logStatistics();
  }

  @Override
  public int numberWindows() {
    return 0;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mReadLength >= 1 && mReadLength <= LONG_BITS);
    Exam.assertTrue(mReadLengthReverse >= 0 && mReadLengthReverse < LONG_BITS);
    Exam.assertTrue("windowSize=" + mWindowSize + " readLength=" + mReadLength, mWindowSize >= 1 && mWindowSize <= mReadLength && mWindowSize <= INT_BITS);
    Exam.assertTrue(mWindowSize * 2 == mWindowBits);
    Exam.assertTrue(bitFlip(mValuesF0) == mValuesR0);
    Exam.assertTrue(bitFlip(mValuesF1) == mValuesR1);
    Exam.assertTrue(mSoFar >= 0);
    if (mTemplateCall instanceof TemplateCallImplementation) {
      final NgsHashFunction hf = ((TemplateCallImplementation) mTemplateCall).mHashFunction;
      Exam.assertTrue(hf == this);
      Exam.integrity(mTemplateCall);
    }
    return true;
  }
}

