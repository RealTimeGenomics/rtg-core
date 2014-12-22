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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.format.FormatReal;

/**
 * Provides some support for the two timer classes
 * <code>(Timer, OneShotTimer)</code> to ensure that they use the
 * same output formats.
 *
 */
@TestClass(value = {"com.rtg.util.diagnostic.TimerTest"})
public abstract class AbstractTimer {

  /** Used as start of output to log (enables timing information to be easily extracted from logs). */
  protected static final String TIMER_PREFIX = "Timer ";

  protected static final FormatReal TIME_FORMAT = new FormatReal(7, 2);

  /** Divisor for calculating seconds from nano-seconds. */
  protected static final double BILLION = 1000000000.0;

  /** The name of the timer (not permitted to contain any spaces). */
  protected String mName;

  /**
   * Create a new timer and start ticking.
   * @param name of times is used in printing timing results.
   */
  public AbstractTimer(final String name) {
    if (name != null && name.contains(" ")) {
      throw new IllegalArgumentException("Name contains spaces:" + name);
    }
    mName = name;
  }

  /**
   * Set name of timer.
   * @param name of timer.
   */
  public void reset(final String name) {
    mName = name;
  }

}

