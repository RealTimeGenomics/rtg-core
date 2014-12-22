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
package com.rtg.vcf.header;

import java.util.HashMap;

/**
 * Encapsulate type in <code>VCF</code> meta lines
 */
public enum MetaType {
  /** if value is integer */
  INTEGER("Integer"),
  /** if value is floating point number */
  FLOAT("Float"),
  /** if value is not present (i.e. existence of ID is sufficient) */
  FLAG("Flag"),
  /** if value is a single character */
  CHARACTER("Character"),
  /** if value is a string */
  STRING("String");

  private static final HashMap<String, MetaType> PARSE_MAP = new HashMap<>(5);
  static {
    for (final MetaType mt : MetaType.values()) {
      PARSE_MAP.put(mt.toString(), mt);
    }
  }
  private final String mToString;

  private MetaType(String toString) {
    mToString = toString;
  }

  /**
   * @return value as appears in file
   */
  @Override
  public String toString() {
    return mToString;
  }

  /**
   * @param val value as appears in file
   * @return corresponding instance
   */
  public static MetaType parseValue(String val) {
    return PARSE_MAP.get(val);
  }

  /**
   * find out if type is super set of another. Currently each type is a super set of itself, and {@link MetaType#FLOAT} is also a super set of {@link MetaType#INTEGER}
   * @param other type to compare to
   * @return true if this is a super set
   */
  public boolean isSuperSet(MetaType other) {
    if (this == MetaType.FLOAT) {
      return other == MetaType.FLOAT || other == MetaType.INTEGER;
    }
    return this == other;
  }
}
