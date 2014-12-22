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
