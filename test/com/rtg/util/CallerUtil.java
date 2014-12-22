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
 * Utitlities for detecting who is calling. Needed in some testing code to
 * avoid bugs in the calling software.
 */
public final class CallerUtil {

  private CallerUtil() { } //so cant instantiate.

  /**
   * Check if being called by the <code>emma</code> coverage tool.
   * @return true iff being called by the <code>emma</code> coverage tool.
   */
  public static boolean emma() {
    return caller("com.vladium.emma");
  }

  /**
   * Check if being called from software whose package begins with name.
   * @param name to be checked.
   * @return true iff name occurs as a prefix on the call stack.
   */
  private static boolean caller(final String name) {
    //
    final StackTraceElement[] trace;
    trace = new Throwable().getStackTrace();
    for (StackTraceElement aTrace : trace) {
      //System.err.println(trace[i].getClassName());
      if (aTrace.getClassName().startsWith(name)) {
        return true;
      }
    }
    return false;
  }


}

