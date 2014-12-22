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

/**
 * Mocks everything except the length of the sequence.
 */
public class ReaderLongMock implements SequencesReader {
  private final long mLength;
  /**
   * @param length length.
   */
  public ReaderLongMock(final long length) {
    mLength = length;
  }
  @Override
  public void close() {
  }
  @Override
  public int currentLength() {
    return 0;
  }
  @Override
  public String currentName() {
    return null;
  }
  @Override
  public String currentFullName() {
    return currentName();
  }
  @Override
  public long currentSequenceId() {
    return 0;
  }
  @Override
  public File path() {
    return null;
  }
  @Override
  public long maxLength() {
    return mLength;
  }
  @Override
  public long minLength() {
    return 0;
  }
  @Override
  public boolean nextSequence() {
    return false;
  }
  @Override
  public long numberSequences() {
    return 1;
  }
  @Override
  public int read(final long sequenceIndex, final byte[] dataOut) {
    return 0;
  }
  @Override
  public int read(long index, byte[] out, int start, int length) {
    return 0;
  }
  @Override
  public int length(long index) {
    return 0;
  }
  @Override
  public String name(long index) {
    return null;
  }
  @Override
  public String fullName(long index) {
    return name(index);
  }
  @Override
  public int readCurrent(byte[] out, int start, int length) {
    return 0;
  }
  @Override
  public int readCurrent(final byte[] dataOut) throws IllegalArgumentException, IllegalStateException {
    return 0;
  }
  @Override
  public void seek(final long sequenceId) {
  }
  @Override
  public long totalLength() {
    return 0;
  }
  @Override
  public SequenceType type() {
    return null;
  }
  @Override
  public long[] residueCounts() {
    return null;
  }
  @Override
  public long dataChecksum() {
    return 0;
  }
  @Override
  public long qualityChecksum() {
    return 0;
  }
  @Override
  public long nameChecksum() {
    return 0;
  }

  private static class NullNames extends PrereadNames {
    @Override
    public String name(final long id) {
      return null;
    }
  }

  @Override
  public PrereadNames names() {
    return new NullNames();
  }

  @Override
  public long lengthBetween(final long start, final long end) {
    return mLength;
  }

  @Override
  public int[] sequenceLengths(final long start, final long end) {
    return new int[] {(int) mLength};
  }

  @Override
  public boolean hasQualityData() {
    return false;
  }

  @Override
  public boolean hasNames() {
    return true;
  }

  @Override
  public int readCurrentQuality(final byte[] dest) throws IllegalArgumentException, IllegalStateException {
    return 0;
  }

  @Override
  public int readCurrentQuality(final byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException {
    return 0;
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest) {
    return 0;
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) {
    return 0;
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

  @Override
  public SdfId getSdfId() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public double[] positionQualityAverage() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long sdfVersion() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public SequencesReader copy() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean compressed() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public void reset() {

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
