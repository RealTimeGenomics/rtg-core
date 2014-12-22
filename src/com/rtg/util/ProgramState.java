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

import com.rtg.util.diagnostic.Diagnostic;

/**
 * <p>Contains static information about program state.
 *
 * <p>Code paths using these mechanisms should be self contained i.e. {@link #setAbort()}
 * and {@link #clearAbort()} should only be called by the main thread and any potential {@link #setAbort()}
 * call should have a <code>finally</code> condition specifying {@link #clearAbort()}
 * once the other threads have safely exited. Note there will be issues for nested structures
 * in this setup.
 *
 */
public final class ProgramState {

  /**
   * Exception thrown when aborting
   */
  public static class SlimAbortException extends RuntimeException {
    /**
     * Constructs the exception
     * @param message error message
     */
    public SlimAbortException(String message) {
      super(message);
    }
  }

  private ProgramState() {

  }

  private static volatile boolean sAbort = false;

  /**
   * Calls to {@link #checkAbort()} after this will throw exception
   */
  public static void setAbort() {
    sAbort = true;
  }

  /**
   * Calls to {@link #checkAbort()} after this will not throw exception
   */
  public static void clearAbort() {
    sAbort = false;
  }

  /**
   * @return true if the abort status has been set by the {@link #setAbort()} method. This method does not throw
   * a {@link SlimAbortException}.
   */
  public static boolean isAbort() {
    return sAbort;
  }

  /**
   * If abort status has been set by the {@link #setAbort()} method then we will throw
   * a {@link SlimAbortException} otherwise we return normally
   */
  public static void checkAbort() {
    if (sAbort) {
      final String message = "Aborting operation in thread: " + Thread.currentThread().getName();
      Diagnostic.userLog(message);
      throw new SlimAbortException(message);
    }
  }
}
