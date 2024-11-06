/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.util.array;


import com.rtg.util.array.bitindex.BitIndex;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;

/**
 * Time doing a "sequential" write of the various IntIndex implementations and some other equivalent methods.
 */
public final class ArrayTiming {

  /** Number of iterations. */
  public static final int ITERATIONS = 1000 * 1000 * 1000;

  /**
   * Get an array type that is capable of storing the specified number of bits
   * in each value. It tries to return the smallest possible result.
   * @param bits to be stored in each unit.
   * @param length number of items to be stored.
   * @param chunked the type of array to use inside the index
   * @return an array type capable of storing the requested number of bits.
   */
  public static CommonIndex bestForBits(final int bits, final long length, final BitIndex.IndexType chunked) {
    if (bits > 64) {
      throw new RuntimeException("Cannot store more than 64 bits per entry bits=" + bits);
    }
    if (bits > 32) {
      if (chunked == BitIndex.IndexType.CHUNKED) {
        return new com.rtg.util.array.longindex.LongChunks(length);
      } else {
        return new com.rtg.util.array.longindex.LongArray(length);
      }
    }
    if (bits > 21) {
      if (chunked == BitIndex.IndexType.CHUNKED) {
        return new com.rtg.util.array.intindex.IntChunks(length);
      } else {
        return new com.rtg.util.array.intindex.IntArray(length);
      }
    }
    if (bits > 16) {
      return new com.rtg.util.array.packedindex.PackedIndex(length, 1L << bits, chunked);
    }
    if (bits == 8) {
      if (chunked == BitIndex.IndexType.CHUNKED) {
        return new com.rtg.util.array.byteindex.ByteChunks(length);
      } else {
        return new com.rtg.util.array.byteindex.ByteArray(length);
      }
    }
    if (bits > 12) {
      if (chunked == BitIndex.IndexType.CHUNKED) {
        return new com.rtg.util.array.shortindex.ShortChunks(length);
      } else {
        return new com.rtg.util.array.shortindex.ShortArray(length);
      }
    }
    if (bits >= 1) {
      return new com.rtg.util.array.packedindex.PackedIndex(length, 1L << bits, chunked);
    }
    if (bits == 0) {
      return new com.rtg.util.array.packedindex.PackedIndex(length, 2L, chunked);
    }
    throw new RuntimeException("bits should be positive:" + bits);
  }

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

  private ArrayTiming() { }

  /**
   * A timing message.
   *  @param msg message
   * @param times array containing experiment times
   */
  public static void tMsg(final String msg, final long[] times) {
    final int lent = times.length;
    tM(times[lent - 1] - times[0]);
    for (int i = 1; i < times.length; ++i) {
      tM(times[i] - times[i - 1]);
    }
    System.out.println("  " + msg);
  }

  private static void tM(final long d) {
    System.out.print(TIME_FORMAT.format(d * 1000000.0 / ITERATIONS));
  }

  /**
   * Time a null loop
   * @param tm time prior to creation
   * @param len IntIndex to be timed
   */
  public static void timeNull(final long tm, final long len) {
    final long t0 = System.currentTimeMillis();
    long j = 0;
    for (int i = 0; i < ITERATIONS; ++i) {
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
   * Time a loop of get/set calls over the index.
   * @param tm time prior to creation
   * @param a IntIndex to be timed
   */
  private static void time(final long tm, final CommonIndex a) {
    if (a != null) {
      final long len = a.length();
      final long t0 = System.currentTimeMillis();
      final long step = 1039111 % len;
      long j = 0;
      for (int i = 0; i < ITERATIONS; ++i) {
        j += step;
        if (j >= len) {
          j = j - len;
        }
        a.get(j);
        a.set(j, 11);
        a.get(j);
      }
      final long t1 = System.currentTimeMillis();
      //a.close();
      final long t2 = System.currentTimeMillis();
      System.gc();
      final long t3 = System.currentTimeMillis();
      tMsg(a.getClass().getName(), new long[] {tm, t0, t1, t2, t3});
    }
  }

  /**
   * Time a loop of get/set calls over the index.
   * @param tm time prior to creation
   * @param a IntIndex to be timed
   * @param threads number of threads
   */
  private static void time(final long tm, final CommonIndex a, final int threads) {
    if (a != null) {
      final Thread[] t = new Thread[threads];
      final long t0 = System.currentTimeMillis();
      for (int k = 0; k < threads; ++k) {
        final int j = k;
        t[k] = new Thread() {

          final long mStartPoint = j;

          @Override
          public void run() {
            final long len = a.length();
            final long step = 1039111 % len;
            long j = mStartPoint;
            for (int i = 0; i < ITERATIONS; ++i) {
              j += step;
              if (j >= len) {
                j = j - len;
              }
              a.get(j);
              a.set(j, 11);
              a.get(j);
            }
          }
        };
      }
      for (final Thread th : t) {
        th.start();
      }
      try {
        for (final Thread th : t) {
          th.join();
        }
      } catch (final InterruptedException e) {
        System.out.println("Interrupted");
      }
      final long t1 = System.currentTimeMillis();
      final long t2 = System.currentTimeMillis();
      System.gc();
      final long t3 = System.currentTimeMillis();
      tMsg(a.getClass().getName(), new long[] {tm, t0, t1, t2, t3});
    }
  }

  /**
   * Print the headings.
   * @param len length
   */
  public static void printHeadings(final long len) {
    System.out.println(LENGTH_FORMAT.format(len) + ":Length  Timing (nanosecs)");
    for (final String element : HEADINGS) {
      System.out.print(element);
    }
    System.out.println();
  }

  /**
   * Run timing tests for the various IntIndex implementations.
   * @param args length (test array) bits threads.
   */
  public static void main(final String[] args) {
    //System.gc();
    if (args.length != 3) {
      System.out.println("Usage: ArrayTiming length bits threads");
      System.exit(1);
    }
    final long len = Long.parseLong(args[0]);
    final int bits  = Integer.parseInt(args[1]);
    final int threads = Integer.parseInt(args[2]);

    printHeadings(len);

    //do timing tests
    timeNull(System.currentTimeMillis(), len);
    timeNull(System.currentTimeMillis(), len);
    timeNull(System.currentTimeMillis(), len);
    try {
      System.out.println("Chunked...");
      for (int i = 0; i < 3; ++i) {
        final CommonIndex index0 = bestForBits(bits, len, BitIndex.IndexType.CHUNKED);
        timeIt(threads, index0);
      }
      System.out.println("Single...");
      for (int i = 0; i < 3; ++i) {
        final CommonIndex index1 = bestForBits(bits, len, BitIndex.IndexType.SINGLE);
        timeIt(threads, index1);
      }
    } catch (final Throwable e) { //might get out of memory
      System.out.println("Exception while timing. " + e + ":" + e.getMessage());
      throw new RuntimeException(e);
    }
  }

  private static void timeIt(final int threads, final CommonIndex index) {
    final long tm = System.currentTimeMillis();
    if (threads > 1) {
      time(tm, index, threads);
    } else {
      time(tm, index);
    }
  }
}
