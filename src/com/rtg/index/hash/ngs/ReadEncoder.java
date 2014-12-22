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

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 * A ReadEncoder encodes read id's during ngs build phase.
 */
public abstract class ReadEncoder implements PseudoEnum {
  /**
   * Paired-end first side.
   */
  public static final ReadEncoder PAIRED_FIRST = new ReadEncoder(0, "PAIRED_FIRST") {
    @Override
    public int encode(final int readId) {
      return readId << 1;
    }
  };

  /**
   * Paired-end second side.
   */
  public static final ReadEncoder PAIRED_SECOND = new ReadEncoder(1, "PAIRED_SECOND") {
    @Override
    public int encode(final int readId) {
      return (readId << 1) + 1;
    }
  };
  /**
   * Single-end reads.
   */
  public static final ReadEncoder SINGLE_END = new ReadEncoder(2, "SINGLE_END") {
    @Override
    public int encode(final int readId) {
      return readId;
    }
  };

  private String mName;
  private int mOrdinal;

  private ReadEncoder(final int ordinal, final String name) {
    mName = name;
    mOrdinal = ordinal;
  }

  /**
   * Encode the read id
   * @param readId the read id
   * @return the encoded value
   */
  public abstract int encode(int readId);

  @Override
  public String toString() {
    return mName;
  }

  @Override
  public String name() {
    return mName;
  }
  @Override
  public int ordinal() {
    return mOrdinal;
  }

  private static final EnumHelper<ReadEncoder> HELPER = new EnumHelper<>(ReadEncoder.class, new ReadEncoder[] {PAIRED_FIRST, PAIRED_SECOND, SINGLE_END});

  /**
   * @return list of the enum values
   */
  public static ReadEncoder[] values() {
    return HELPER.values();
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str the name of the enum
   * @return the enum value
   */
  public static ReadEncoder valueOf(final String str) {
    return HELPER.valueOf(str);
  }
}

