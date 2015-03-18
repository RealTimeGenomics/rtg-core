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
package com.rtg.reader;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.mode.SequenceType;

/**
 * Dummy needed for constructing tests.
 */
public class MockSequencesReader extends DummySequencesReader {

  private static final class SimplePrereadNames extends PrereadNames {
    @Override
    public String name(final long id) {
      return "seq" + id;
    }
  }

  private final SequenceType mSequenceType;

  private final long mNumberSequences;

  protected long mTotalLength;

  private int[] mLengths;

  /**
   * @param sequenceType sequence type.
   * @param numberSquences number sequences.
   */
  public MockSequencesReader(final SequenceType sequenceType, final long numberSquences) {
    mSequenceType = sequenceType;
    mNumberSequences = numberSquences;
    mTotalLength = 0;
    if (mNumberSequences <= Integer.MAX_VALUE) {
      mLengths = new int[(int) mNumberSequences];
      Arrays.fill(mLengths, 1);
    }
  }

  /**
   * @param sequenceType sequence type.
   * @param numberSquences number sequences.
   * @param length total length of the sequences.
   */
  public MockSequencesReader(final SequenceType sequenceType, final long numberSquences, final long length) {
    mSequenceType = sequenceType;
    mNumberSequences = numberSquences;
    mTotalLength = length;
    if (mNumberSequences <= Integer.MAX_VALUE) {
      mLengths = new int[(int) mNumberSequences];
      if (mNumberSequences > 0) {
        final int fillVal = (int) (mTotalLength / mNumberSequences);
        Arrays.fill(mLengths, fillVal);
        mLengths[0] += mTotalLength % mNumberSequences;
      }
    }
  }

  /**
   * @param sequenceType sequence type.
   */
  public MockSequencesReader(final SequenceType sequenceType) {
    mSequenceType = sequenceType;
    mNumberSequences = 100;
    mTotalLength = 0;
    mLengths = new int[(int) mNumberSequences];
    Arrays.fill(mLengths, 1);
  }

  /**
   * Set an array of lengths for returning
   * @param lengths the lengths of the pseudo sequences
   */
  public void setLengths(int[] lengths) {
    if (lengths != null) {
      mLengths = Arrays.copyOf(lengths, lengths.length);
      mTotalLength = 0;
      for (final int length : mLengths) {
        mTotalLength += length;
      }
    } else {
      mLengths = null;
    }
  }

  @Override
  public long totalLength() {
    return mTotalLength;
  }

  @Override
  public long maxLength() {
    return 1000;
  }
  @Override
  public long minLength() {
    return 1;
  }
  @Override
  public long numberSequences() {
    return mNumberSequences;
  }

  @Override
  public SequenceType type() {
    return mSequenceType;
  }

  @Override
  public int length(long index) {
    if (mLengths != null && index >= 0 && index < mLengths.length) {
      return mLengths[(int) index];
    }
    return 0;
  }
  @Override
  public String name(long index) {
    return "seq" + index;
  }
  @Override
  public String fullName(long index) {
    return name(index);
  }

  @Override
  public PrereadNames names() {
    return new SimplePrereadNames();
  }

  @Override
  public long lengthBetween(final long start, final long end) {
    long tot = 0;
    for (long i = start; i < end; i++) {
      tot += mLengths[(int) i];
    }
    return tot;
  }

  @Override
  public int[] sequenceLengths(final long start, final long end) {
    final int[] ret = new int[(int) (end - start)];
    System.arraycopy(mLengths, (int) start, ret, 0, (int) (end - start));
    return ret;
  }

  @Override
  public boolean hasNames() {
    return true;
  }

  @Override
  public int readQuality(final long sequenceIndex, final byte[] dest) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean hasHistogram() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long[] histogram() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long[] posHistogram() {
    throw new UnsupportedOperationException("Not in memory readers");
  }

  @Override
  public long longestNBlock() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long nBlockCount() {
    throw new UnsupportedOperationException("Not supported yet.");
  }


  @Override
  public double globalQualityAverage() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public PrereadArm getArm() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public PrereadType getPrereadType() {
    throw new UnsupportedOperationException("Not supported yet.");
  }


}
