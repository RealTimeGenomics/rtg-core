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
package com.rtg.position.output;

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 * Enum class
 */
public final class Offset implements PseudoEnum {

  /** Enum member */
  public static final Offset ZERO = new Offset(0, "ZERO");
  /** Enum member */
  public static final Offset INSERT = new Offset(1, "INSERT");
  /** Enum member */
  public static final Offset DELETE = new Offset(2, "DELETE");
  /** Enum member */
  public static final Offset TOTAL = new Offset(3, "TOTAL");
  private final String mName;
  private final int mOrdinal;

  private Offset(final int ordinal, final String name) {
    super();
    mOrdinal = ordinal;
    mName = name;
  }

  @Override
  public int ordinal() {
    return mOrdinal;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public String toString() {
    return mName;
  }

  private static final EnumHelper<Offset> HELPER = new EnumHelper<>(Offset.class, new Offset[]{ZERO, INSERT, DELETE, TOTAL});

  /**
   * {@link EnumHelper#valueOf(String)}
   * @param str as in EnumHelper
   * @return as in EnumHelper
   */
  public static Offset valueOf(final String str) {
    return HELPER.valueOf(str);
  }

  /**
   * {@link EnumHelper#values()}
   * @return as in EnumHelper
   */
  public static Offset[] values() {
    return HELPER.values();
  }
}
