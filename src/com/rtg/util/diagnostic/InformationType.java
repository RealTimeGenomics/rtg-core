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
 * Enumeration of SLIM information messages.
 * See <code>src.com.reeltwo.cartesian.util.diagnostic.Diagnostics.properties</code>
 * for the localised messages.
 *
 */
public final class InformationType implements DiagnosticType, PseudoEnum, Serializable {

  private static int sCounter = -1;

  /**
   * Simple User Message
   * You should supply one argument which is the message to display
   */
  public static final InformationType INFO_USER = new InformationType(++sCounter, "INFO_USER", 1);

  /**
   * Information about what file is being processed
   * <code>'Processing %1"%2" (%3 of %4)'</code>
   */
  public static final InformationType PROCESSING_ITEM_N_OF_N = new InformationType(++sCounter, "PROCESSING_ITEM_N_OF_N", 4);

  private static final EnumHelper<InformationType> HELPER = new EnumHelper<>(InformationType.class, new InformationType[] {
    INFO_USER,
    PROCESSING_ITEM_N_OF_N
  });

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str name of value
   * @return the enum value
   */
  public static InformationType valueOf(final String str) {
    return HELPER.valueOf(str);
  }

  /**
   * @return list of enum values
   */
  public static InformationType[] values() {
    return HELPER.values();
  }


  /** Number of parameters that must occur in conjunction with this information. */
  private final int mParams;

  private final int mOrdinal;

  private final String mName;

  private InformationType(final int ordinal, final String name, final int params) {
    mParams = params;
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

  @Override
  public int getNumberOfParameters() {
    return mParams;
  }

  Object readResolve() {
    return values()[this.ordinal()];
  }

  @Override
  public String getMessagePrefix() {
    return "";
  }

}
