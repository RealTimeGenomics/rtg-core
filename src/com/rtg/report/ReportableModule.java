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

package com.rtg.report;

/**
 */
public enum ReportableModule {
  /** species module */
  SPECIES(10),
  /** similarity module */
  SIMILARITY(1),
  /** map module */
  MAP(500);

  private final int mMaxInputDirectories;

  /**
   * Constructor for a maximum number of input directories for a report type.
   * @param maxInputDirectories maximum input directories for this report type.
   */
  private ReportableModule(int maxInputDirectories) {
    mMaxInputDirectories = maxInputDirectories;
  }

  /**
   * Get the maximum number of input directories allowed for this report type.
   * @return the maximum number of input directories allowed for this report type.
   */
  public int getMaxInputDirectories() {
    return mMaxInputDirectories;
  }
}
