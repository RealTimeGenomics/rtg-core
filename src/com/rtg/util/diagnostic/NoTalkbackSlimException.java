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

/**
 * Exception for reporting errors within SLIM which are a result of local user
 * or machine error.  The exceptions are logged but are not reported by the talkback mechanism.
 * Being a subclass of
 * {@link SlimException} they allow the program to terminate with the
 * proper error code.
 *
 */
public class NoTalkbackSlimException extends SlimException {

  /**
   * A new exception for an error where intelligent information can be
   * returned to the user via the diagnostics system.
   * A stack trace will not be logged.
   *
   * @param type error type
   * @param args parameters of the error type
   */
  public NoTalkbackSlimException(final ErrorType type, final String... args) {
    this(null, type, args);
  }

  /**
   * A new exception for an error where intelligent information can be
   * returned to the user via the diagnostics system.
   * A stack trace will not be logged.
   *
   * @param t ignored
   * @param type error type
   * @param args parameters of the error type
   */
  public NoTalkbackSlimException(final Throwable t, final ErrorType type, final String... args) {
    super(false, t, type, args);
  }

  /**
   * Report program failure to user
   * @param message message for error
   */
  public NoTalkbackSlimException(String message) {
    this(null, ErrorType.INFO_ERROR, message);
  }
}
