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
 * Provide output to Diagnostic of timer information.
 * This timer is very simple. All it can do is be created and
 * log the time. It uses the same conventions for output format as
 * <code>Timer</code>.
 *
 * Precision is to nano-seconds and output is in fractions of a second.
 *
 * Output to the log is done in a standard form with the word "Time" as the
 * first token followed by an identifier (without any embedded spaces)
 * followed by the total time.
 *
 */
public class OneShotTimer extends AbstractTimer {

  /** The time at the start. */
  private final long mStart;

  /**
   * Create a new timer and start ticking.
   * @param name name of timer
   */
  public OneShotTimer(final String name) {
    super(name);
    mStart = System.nanoTime();
  }

  /**
   * Finish timing and log the resulting interval.
   */
  public void stopLog() {
    final long finish = System.nanoTime();
    final long diff = finish - mStart;
    Diagnostic.developerLog(toString(diff));
  }

  /**
   * @param diff time difference in nanoseconds to be displayed.
   * @return formatted time output using time difference.
   */
  String toString(final long diff) {
    final StringBuilder sb = new StringBuilder();
    sb.append(TIMER_PREFIX).append(mName);
    //System.err.println("start=" + mStart + " finish=" + finish);
    final double time = diff / BILLION;
    TIME_FORMAT.format(sb, time);
    return sb.toString();
  }

  @Override
  public String toString() {
    return mName;
  }
}

