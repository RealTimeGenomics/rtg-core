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

import com.rtg.mode.DnaUtils;

/**
 */
public class ArraySequencesReader extends DummySequencesReader {
  private final byte[][] mData;
  private final byte[][] mQuality;
  private final int mTotalLength;

  /**
   * @param data strings containing ascii DNA sequence
   */
  public ArraySequencesReader(final String... data) {
    this(convertStrings(data), null);
  }

  private static byte[][] convertStrings(String... data) {
    final byte[][] result = new byte[data.length][];
    for (int i = 0; i < data.length; i++) {
      result[i] = DnaUtils.encodeString(data[i]);
    }
    return result;
  }

  /**
   * @param data sequences of nucleotides (0 to 4 convention).
   * @param quality quality for sequences (can be null).
   */
  public ArraySequencesReader(final byte[][] data, final byte[][] quality) {
    mData = new byte[data.length][];
    mQuality = quality == null ? null : new byte[quality.length][];
    int tot = 0;
    for (int i = 0; i < data.length; i++) {
      final byte[] element = data[i];
      mData[i] = element.clone();
      tot += element.length;
    }
    if (mQuality != null) {
      for (int i = 0; i < quality.length; i++) {
        mQuality[i] = quality[i].clone();
      }
    }
    mTotalLength = tot;
  }

  @Override
  public long totalLength() {
    return mTotalLength;
  }

  @Override
  public long maxLength() {
    long max = 0;
    for (final byte[] element : mData) {
      if (element.length > max) {
        max = element.length;
      }
    }
    return max;
  }

  @Override
  public long minLength() {
    long min = Integer.MAX_VALUE;
    for (final byte[] element : mData) {
      if (element.length < min) {
        min = element.length;
      }
    }
    return min;
  }

  @Override
  public long numberSequences() {
    return mData.length;
  }


  @Override
  public byte[] read(final long sequenceIndex) {
    final byte[] result = new byte[length(sequenceIndex)];
    read(sequenceIndex, result);
    return result;
  }

  @Override
  public int read(final long sequenceIndex, final byte[] dataOut) {
    final byte[] data = mData[(int) sequenceIndex];
    final int length = data.length;
    if (length > dataOut.length) {
      throw new IllegalArgumentException();
    }
    System.arraycopy(data, 0, dataOut, 0, length);
    return length;
  }

  @Override
  public int read(final long index, final byte[] out, final int start, final int length) {
    final byte[] data = mData[(int) index];
    final int lengthData = data.length;
    if (start + length > lengthData || length > out.length) {
      throw new IllegalArgumentException();
    }
    System.arraycopy(data, start, out, 0, length);
    return length;
  }

  @Override
  public String name(long index) {
    return "sequence " + index;
  }

  @Override
  public String fullName(long index) {
    return name(index);
  }

  @Override
  public int length(final long index) {
    return mData[(int) index].length;
  }

  @Override
  public long lengthBetween(final long start, final long end) {
    long tot = 0;
    for (int i = (int) start; i < end; i++) {
      tot += mData[i].length;
    }
    return tot;
  }

  @Override
  public int[] sequenceLengths(final long start, final long end) {
    final int[] lengths = new int[mData.length];
    for (int i = 0; i < mData.length; i++) {
      lengths[i] = mData[i].length;
    }
    return lengths;
  }

  @Override
  public boolean hasQualityData() {
    return mQuality != null;
  }

  @Override
  public boolean hasNames() {
    return true;
  }

  @Override
  public byte[] readQuality(final long sequenceIndex) {
    final byte[] result = new byte[mQuality[(int) sequenceIndex].length];
    readQuality(sequenceIndex, result);
    return result;
  }

  @Override
  public int readQuality(final long sequenceIndex, final byte[] dest) {
    return readQuality(sequenceIndex, dest, 0, mQuality[(int) sequenceIndex].length);
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) {
    if (mQuality == null) {
      throw new IllegalStateException();
    }
    final byte[] quality = mQuality[(int) sequenceIndex];
    if (start + length > dest.length) {
      throw new IllegalArgumentException();
    }
    System.arraycopy(quality, start, dest, 0, length);
    return length;
  }

}
