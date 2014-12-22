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
package com.rtg.launcher;


import java.util.Locale;
import java.util.ResourceBundle;

/**
 * Common static information and methods for the various build programs.
 */
public final class BuildCommon {

  private BuildCommon() { }

  /** Properties file with all non-logged strings used in the system. */
  static final String BUILDCOMMON_RESOURCE_BUNDLE = "com.rtg.launcher.BuildCommon";
  /** Internationalized messages. */
  public static final ResourceBundle RESOURCE = ResourceBundle.getBundle(BUILDCOMMON_RESOURCE_BUNDLE, Locale.getDefault());



  /** Query flag. */
  public static final String QUERY_FLAG = RESOURCE.getString("QUERY_FLAG");
  /** Subject flag. */
  public static final String SUBJECT_FLAG = RESOURCE.getString("SUBJECT_FLAG");
  /** Program mode flag. */
  public static final String PROGRAM_MODE_FLAG = RESOURCE.getString("PROGRAM_MODE_FLAG");
  /** Maximum allowed gap. */
  public static final String MAXGAP_FLAG = RESOURCE.getString("MAXGAP_FLAG");
  /** Progress flag. */
  public static final String PROGRESS_FLAG = RESOURCE.getString("PROGRESS_FLAG");

}

