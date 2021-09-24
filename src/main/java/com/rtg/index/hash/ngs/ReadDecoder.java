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
 *
 */
public abstract class ReadDecoder implements PseudoEnum {
  /**
   * Paired-end mode.
   */
  public static final ReadDecoder PAIRED_END = new ReadDecoder(0, "PAIRED_END") {
    @Override
    public int decode(final int internalId) {
      return internalId >>> 1;
    }

    @Override
    public boolean isFirst(int internalId) {
      return (internalId & 1) == 0;
    }
  };

  /**
   * Single-end reads.
   */
  public static final ReadDecoder SINGLE_END = new ReadDecoder(1, "SINGLE_END") {
    @Override
    public int decode(final int internalId) {
      return internalId;
    }

    @Override
    public boolean isFirst(int internalId) {
      return true;
    }
  };

  private final String mName;
  private final int mOrdinal;

  private ReadDecoder(final int ordinal, final String name) {
    mName = name;
    mOrdinal = ordinal;
  }

  /**
   * Decode the internal id to a read id
   * @param internalId the internal id
   * @return the decoded read id
   */
  public abstract int decode(int internalId);

  /**
   * Whether this internal id is from the first or second arm
   * @param internalId the internal id
   * @return true if first
   */
  public abstract boolean isFirst(int internalId);
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

  private static final EnumHelper<ReadDecoder> HELPER = new EnumHelper<>(ReadDecoder.class, new ReadDecoder[] {PAIRED_END, SINGLE_END});

  /**
   * @return list of the enum values
   */
  public static ReadDecoder[] values() {
    return HELPER.values();
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str the name of the enum
   * @return the enum value
   */
  public static ReadDecoder valueOf(final String str) {
    return HELPER.valueOf(str);
  }
}
