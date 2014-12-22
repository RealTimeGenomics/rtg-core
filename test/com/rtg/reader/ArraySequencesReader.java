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

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceType;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class ArraySequencesReader extends IntegralAbstract implements SequencesReader {
  private final byte[][] mData;
  private final byte[][] mQuality;
  private final int mTotalLength;

  private int mIndex;

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
    mIndex = -1;
  }

  @Override
  public void seek(final long sequenceId) {
    mIndex = (int) sequenceId;
    assert integrity();
  }

  @Override
  public boolean nextSequence() {
    mIndex++;
    return mIndex < mData.length;
  }

  @Override
  public SequenceType type() {
    throw new UnsupportedOperationException("Not implemented yet.");
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
  public long currentSequenceId() {
    return mIndex;
  }

  @Override
  public long dataChecksum() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }
  @Override
  public long qualityChecksum() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }
  @Override
  public long nameChecksum() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public long[] residueCounts() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public long[] histogram() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public long[] posHistogram() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public double globalQualityAverage() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public double[] positionQualityAverage() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public long nBlockCount() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public long longestNBlock() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public boolean hasHistogram() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public PrereadArm getArm() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public PrereadType getPrereadType() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public SdfId getSdfId() {
    throw new UnsupportedOperationException("Not implemented yet.");
  }

  @Override
  public int currentLength() {
    return mData[mIndex].length;
  }

  @Override
  public String currentName() {
    return "sequence " + mIndex;
  }

  @Override
  public String currentFullName() {
    return currentName();
  }

  @Override
  public SequencesReader copy() {
    throw new UnsupportedOperationException("Not supported yet.");
  }


  @Override
  public File path() {
    throw new UnsupportedOperationException("Not implemented yet.");
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
  public int readCurrent(final byte[] out, final int start, final int length) {
    return read(mIndex, out, start, length);
  }

  @Override
  public int readCurrent(final byte[] dataOut) throws IllegalArgumentException, IllegalStateException {
    return read(mIndex, dataOut);
  }

  @Override
  public void close() {
    //do nothing
  }

  @Override
  public PrereadNames names() {
    throw new UnsupportedOperationException("Not implemented yet.");
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
  public int readQuality(final long sequenceIndex, final byte[] dest) {
    return readQuality(sequenceIndex, dest, 0, mQuality[mIndex].length);
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) {
    if (mQuality == null) {
      throw new IllegalStateException();
    }
    final byte[] quality = mQuality[mIndex];
    if (start + length > dest.length) {
      throw new IllegalArgumentException();
    }
    System.arraycopy(quality, start, dest, 0, length);
    return length;
  }

  @Override
  public int readCurrentQuality(final byte[] dest) throws IllegalArgumentException, IllegalStateException {
    return readQuality(mIndex, dest);
  }

  @Override
  public int readCurrentQuality(byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException {
    return readQuality(mIndex, dest, start, length);
  }

  @Override
  public long sdfVersion() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ByteArraySequencesReader");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mData);
    Exam.assertTrue(-1 <= mIndex && mIndex < mData.length);
    Exam.assertTrue(0 <= mTotalLength);
    return true;
  }

  @Override
  public boolean compressed() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void reset() {
    mIndex = -1;
  }

  @Override
  public String currentNameSuffix() throws IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String nameSuffix(long sequenceIndex) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long suffixChecksum() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String getReadMe() {
    return null;
  }

}
