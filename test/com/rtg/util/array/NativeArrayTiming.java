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
package com.rtg.util.array;

import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;

/**
 * Time doing a "sequential" write of the various IntIndex implementations and some other equivalent methods.
 */
public final class NativeArrayTiming {

  /** Number of iterations. */
  public static final int ITERATIONS = 100000000;
  /** Formatting. */
  public static final FormatInteger LENGTH_FORMAT = new FormatInteger(12, true);

  static final FormatReal TIME_FORMAT = new FormatReal(5, 3);

  private static final String[] HEADINGS = {
    "     ",
    "Total    ",
    "Creation ",
    "Loop     ",
    "Close    ",
    "GC       ",
  "type   "};

  private NativeArrayTiming() { }

  /**
   * A timing message.
   *  @param msg message
   * @param times array containing experiment times
   */
  public static void tMsg(final String msg, final long[] times) {
    final int lent = times.length;
    tM(times[lent - 1] - times[0]);
    for (int i = 1; i < times.length; i++) {
      tM(times[i] - times[i - 1]);
    }
    System.out.println("  " + msg);
  }

  private static void tM(final long d) {
    System.out.print(TIME_FORMAT.format(d / 1000.0));
  }

  /**
   * Time a null loop
   * @param tm time prior to creation
   * @param len IntIndex to be timed
   */
  public static void timeNull(final long tm, final long len) {
    final long t0 = System.currentTimeMillis();
    long j = 0;
    for (int i = 0; i < ITERATIONS; i++) {
      j += 1039111L;
      if (j >= len) {
        j = j - len;
      }
    }
    final long t1 = System.currentTimeMillis();
    System.gc();
    final long t3 = System.currentTimeMillis();
    tMsg("null loop", new long[] {tm, t0, t1, t1, t3});
  }

  /**
   * Time a loop of get/set calls over an array.
   * @param tm time prior to creation
   * @param len IntIndex to be timed
   */
  public static void timeArray(final long tm, final int len) {
    final int[] a = new int[len];
    final long t0 = System.currentTimeMillis();
    final int step = (int) (1039111L % len);
    int j = 0;
    for (int i = 0; i < ITERATIONS; i++) {
      j += step;
      if (j >= len) {
        j = j - len;
      }
      a[j] = 11;
    }
    final long t1 = System.currentTimeMillis();
    System.gc();
    final long t3 = System.currentTimeMillis();
    tMsg("array", new long[] {tm, t0, t1, t1, t3});
  }

  /**
   * Print the headings.
   * @param len length
   */
  public static void printHeadings(final long len) {
    System.out.println(LENGTH_FORMAT.format(len) + ":Length  Timing (secs)");
    for (final String element : HEADINGS) {
      System.out.print(element);
    }
    System.out.println();
  }

  /**
   * Run timing tests for the various IntIndex implementations.
   * @param args single argument specifying the length of the indexes to be timed.
   */
  public static void main(final String[] args) {
    //System.gc();
    final long len = Long.parseLong(args[0]) - 1L;
    final int threads = args.length > 2 ? Integer.parseInt(args[2]) : 0;

    printHeadings(len);

    //do timing tests
    timeNull(System.currentTimeMillis(), len);
    timeNull(System.currentTimeMillis(), len);
    timeNull(System.currentTimeMillis(), len);
    if (threads == 1) { //the implementations with a maximum length
      timeArray(System.currentTimeMillis(), (int) len);
      timeArray(System.currentTimeMillis(), (int) len);
      timeArray(System.currentTimeMillis(), (int) len);
    } else {
      try {
        //TODO
      } catch (final Throwable e) { //might get out of memory
        System.out.println("Exception while timing Chunks. " + e.getMessage());
      }
    }
  }
}

