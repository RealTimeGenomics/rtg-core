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
package com.rtg.util;

/**
 * Various constants applying globally to SLIM.
 *
 */
public final class Constants {

  private Constants() { }

  /** Main application name */
  public static final String APPLICATION_NAME = "rtg";

  /**
   * Minimum file chunking size.  This value should be large enough to
   * accommodate a single record in files written by SLIM.
   */
  public static final long MINIMUM_FILE_CHUNK_SIZE = 1000;

  /**
   * Default maximum size that individual files can have
   * (at least the ones we write).
   */
  public static final long MAX_FILE_SIZE = 1000000000L; //1GigaByte

  /** Base for RTG email addresses */
  static final String BASE_EMAIL_ADDR = "@realtimegenomics.com";

  /** Support email address */
  public static final String SUPPORT_EMAIL_ADDR = "support" + BASE_EMAIL_ADDR;

  /** Maximum threads to allow from automatic thread assignment **/
  public static final int MAX_IO_THREADS = 4;

  /** Maximum number of files to open at once **/
  public static final int MAX_OPEN_FILES = 400; //Integer.parseInt(System.getProperty("rtg.max_open_files", "400"));

  /** Number of bytes in a KB */
  public static final double KB = 1024;

  /** Number of bytes in a MB */
  public static final double MB = 1024 * KB;

  /** Number of bytes in a GB */
  public static final double GB = 1024 * MB;
}

