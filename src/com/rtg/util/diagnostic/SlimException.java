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
 * Exception for reporting errors within SLIM.  In every case the exception
 * is logged via the diagnostics system.
 *
 */
public class SlimException extends RuntimeException {

  private final boolean mTalkback;
  private final ErrorType mErrorType;
  private final String[] mErrorArguments;

  /**
   * A new exception for an error where intelligent information can be
   * returned to the user via the diagnostics system and where an
   * underlying exception may contain further information useful for debuggers.
   *
   * @param talkback if true then invoke talkback mechanism.
   * @param t cause
   * @param type error type
   * @param args parameters of the error type
   */
  protected SlimException(final boolean talkback, final Throwable t, final ErrorType type, final String... args) {
    super(t != null ? t.toString() : args.length > 0 ? args[0] : "", t);  //TODO this is heinous and should be dealt to with extreme prejudice
    mTalkback = talkback;
    mErrorType = type;
    mErrorArguments = args.clone();
  }


  /**
   * A new exception for an error where intelligent information can be
   * returned to the user via the diagnostics system and where an
   * underlying exception may contain further information useful for debuggers.
   *
   * @param t cause
   * @param type error type
   * @param args parameters of the error type
   */
  public SlimException(final Throwable t, final ErrorType type, final String... args) {
    this(true, t, type, args);
  }

  /**
   * A new exception for an error where intelligent information can be
   * returned to the user via the diagnostics system.
   * A stack trace will be logged.
   *
   * @param type error type
   * @param args parameters of the error type
   */
  public SlimException(final ErrorType type, final String... args) {
    this(null, type, args);
  }

  /**
   * A new exception wrapping another throwable.
   * A stack trace will be logged.
   *
   * @param t cause
   */
  public SlimException(final Throwable t) {
    this(t, ErrorType.SLIM_ERROR);
  }

  /**
   * A new exception where no additional useful information is available.
   * A stack trace will be logged.
   */
  public SlimException() {
    this(ErrorType.SLIM_ERROR);
  }

  /**
   * An exception where we have an error message for the user.
   * @param message message for the error.
   */
  public SlimException(final String message) {
    this(ErrorType.INFO_ERROR, message);
  }

  /**
   * The error that caused this exception
   * @return the type
   */
  public ErrorType getErrorType() {
    return mErrorType;
  }

  /**
   * Log the error to the diagnostic log file.
   */
  public void logException() {
    Diagnostic.userLog(getCause() == null ? this : getCause());
    Diagnostic.errorLogOnly(mErrorType, mErrorArguments);
  }

  /**
   * Print the user viewable error without logging.
   */
  public void printErrorNoLog() {
    Diagnostic.errorNoLog(mErrorType, mErrorArguments);
  }

  /**
   * Invoke the talkback mechanism if this was a talkback exception.
   */
  public void invokeTalkback() {
    if (mTalkback) {
      Talkback.postTalkback(getCause() == null ? this : getCause());
    }
  }
}

