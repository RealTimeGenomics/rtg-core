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

/**
 * Mocks everything except the length of the sequence.
 */
public class ReaderLongMock extends DummySequencesReader {
  private final long mLength;
  /**
   * @param length length.
   */
  public ReaderLongMock(final long length) {
    mLength = length;
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
  public long numberSequences() {
    return 1;
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
  public boolean hasNames() {
    return true;
  }

}
