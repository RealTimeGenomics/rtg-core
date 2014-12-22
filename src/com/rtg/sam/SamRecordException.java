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
package com.rtg.sam;

/**
 * Indicates problem in the a SAM or BAM record
 */
public class SamRecordException extends RuntimeException {

  /**
   * Constructs exception with given message
   * @param message the message
   */
  public SamRecordException(String message) {
    super(message);
  }

  /**
   * Constructs exception with given message and cause
   * @param message the message
   * @param cause the cause
   */
  public SamRecordException(String message, Throwable cause) {
    super(message, cause);
  }
}
