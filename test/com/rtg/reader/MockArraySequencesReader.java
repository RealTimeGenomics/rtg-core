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

import com.rtg.mode.SequenceType;

/**
 * Mock for constructing tests.
 * The lengths of the reads are specified but will fail if any attempt is made to get at the contents.
 */
public class MockArraySequencesReader extends MockSequencesReader {

  private final String[] mNames;

  private final long mMaxLength;

  private final long mMinLength;


  /**
   * @param sequenceType sequence type.
   * @param length number of sequences.
   */
  public MockArraySequencesReader(final SequenceType sequenceType, final long length) {
    this(sequenceType, length, 1);
  }

  /**
   * @param sequenceType sequence type.
   * @param length number of sequences.
   * @param maxLength maximum length
   */
  public MockArraySequencesReader(final SequenceType sequenceType, final long length, final long maxLength) {
    super(sequenceType, length);
    mNames = null;
    mMaxLength = maxLength;
    mMinLength = -1;
  }

  /**
   * @param sequenceType sequence type.
   * @param seqLengths lengths of sequences.
   */
  public MockArraySequencesReader(final SequenceType sequenceType, final int[] seqLengths) {
    this(sequenceType, seqLengths, null);
  }

  /**
   * @param sequenceType sequence type.
   * @param seqLengths lengths of sequences.
   * @param names names of sequences.
   */
  public MockArraySequencesReader(final SequenceType sequenceType, final int[] seqLengths, final String[] names) {
    super(sequenceType, seqLengths.length, -1);
    mNames = names == null ? null : names.clone();
    super.setLengths(seqLengths);
    long max = 0;
    long min = Integer.MAX_VALUE;
    for (final int length : seqLengths) {
      if (length > max) {
        max = length;
      }
      if (length < min) {
        min = length;
      }
    }
    mMaxLength = max;
    mMinLength = min;
  }

  @Override
  public String name(long index) {
    if (mNames == null) {
      return "seq" + index;
    } else {
      return mNames[(int) index];
    }
  }

  @Override
  public long maxLength() {
    if (mMaxLength == -1) {
      throw new UnsupportedOperationException();
    }
    return mMaxLength;
  }
  @Override
  public long minLength() {
    if (mMinLength == -1) {
      throw new UnsupportedOperationException();
    }
    return mMinLength;
  }

  @Override
  public int read(final long index, final byte[] out) {
    for (int i = 0; i < out.length; i++) {
      out[i] = (byte) ((i % 4) + 1);
    }
    return out.length;
  }

  @Override
  public PrereadNames names() {
    return new PrereadNames() {
        @Override
        public String name(final long id) {
          if (mNames == null) {
            return "seq" + id;
          } else {
            return mNames[(int) id];
          }
        }
      };
  }

}
