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
package com.rtg.util.diagnostic;

import java.io.Serializable;

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 */
public final class ListenerType implements PseudoEnum, Serializable {

  /**
   * Command line listener. Writes messages including progress to stderr.
   */
  public static final ListenerType CLI = new ListenerType(0, "CLI");

  /**
   * File listener. Writes progress messages to the file "progress" in the output directory.
   * The file is rewritten so it contains a single line with the latest progress.
   */
  public static final ListenerType FILE = new ListenerType(1, "FILE");

  /**
   * Specifies that no event listener is to be used.
   */
  public static final ListenerType NULL = new ListenerType(2, "NULL");

  private static final EnumHelper<ListenerType> HELPER = new EnumHelper<>(ListenerType.class, new ListenerType[] {CLI, FILE, NULL});

  /**
   * @return list of the enum values
   */
  public static ListenerType[] values() {
    return HELPER.values();
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str the name of the enum
   * @return the enum value
   */
  public static ListenerType valueOf(final String str) {
    return HELPER.valueOf(str);
  }


  private String mName;
  private int mOrdinal;

  private ListenerType(final int ordinal, final String name) {
    mName = name;
    mOrdinal = ordinal;
  }

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
}
