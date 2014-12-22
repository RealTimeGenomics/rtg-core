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
package com.rtg.util.cli;


/**
 */
public final class CommonFlagCategories {
  private CommonFlagCategories() { }

  /** input output category */
  public static final String INPUT_OUTPUT = "File Input/Output";

  /** sensitivity tuning category */
  public static final String SENSITIVITY_TUNING = "Sensitivity Tuning";

  /** utility category */
  public static final String UTILITY = "Utility";

  /** reporting category */
  public static final String REPORTING = "Reporting";

  /** filtering category */
  public static final String FILTERING = "Filtering";

  /** categories used in display of flags  */
  private static final String[] CATEGORIES = {
    INPUT_OUTPUT,
    SENSITIVITY_TUNING,
    FILTERING,
    REPORTING,
    UTILITY
  };

  /**
   * Sets the categories
   * @param flags flags to set categories on
   */
  public static void setCategories(final CFlags flags) {
    flags.setCategories(UTILITY, CATEGORIES);
  }
}
