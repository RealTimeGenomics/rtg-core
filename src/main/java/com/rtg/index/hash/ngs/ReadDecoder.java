/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
