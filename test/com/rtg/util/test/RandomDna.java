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
package com.rtg.util.test;

import com.rtg.util.PortableRandom;

/**
 * Generate random DNA nucleotides according to the uniform distribution over
 * A, G, C, and T.
 */
public final class RandomDna {

  private RandomDna() { }

  private static final PortableRandom RANDOM = new PortableRandom();

  /**
   * Random nucleotide string.
   *
   * @param len length of string
   * @param r random number generator
   * @return random string
   */
  public static String random(final int len, final PortableRandom r) {
    final StringBuilder sb = new StringBuilder(len);
    for (int i = 0; i < len; i++) {
      switch (r.nextInt(4)) {
      case 0:
        sb.append('A');
        break;
      case 1:
        sb.append('C');
        break;
      case 2:
        sb.append('G');
        break;
      default:
        sb.append('T');
        break;
      }
    }
    return sb.toString();
  }

  /**
   * Random nucleotide string.
   *
   * @param len length of string
   * @return random string
   */
  public static String random(final int len) {
    return random(len, RANDOM);
  }

  /**
   * Main program.
   *
   * @param args size of sequence
   */
  public static void main(final String[] args) {
    if (args.length < 1 && args.length > 2) {
      System.err.println("USAGE: total-size sequence-size");
      return;
    }
    int seqNum = 0;
    if ("-random".equals(args[0])) {
      while (true) {
        System.out.println(">d" + seqNum++);
        System.out.println(random(RANDOM.nextInt(0x000FFFFF)));
      }
    } else {
      long size = Long.parseLong(args[0]);
      long chunk = 1000;
      if (args.length > 1) {
        chunk = Long.parseLong(args[1]);
      }

      while (size > 0) {
        System.out.println(">d" + seqNum++);
        System.out.println(random((int) Math.min(chunk, size)));
        size -= chunk;
      }
    }
  }
}

