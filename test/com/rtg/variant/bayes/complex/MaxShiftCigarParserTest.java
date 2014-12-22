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

package com.rtg.variant.bayes.complex;

import junit.framework.TestCase;

/**
 */
public class MaxShiftCigarParserTest extends TestCase {


  public void testStartPos() {

    final MaxShiftCigarParser msc = new MaxShiftCigarParser();

    msc.parse("1I", 0);
    assertEquals(0, msc.getStartPos());

    msc.parse("2I", 0);
    assertEquals(-1, msc.getStartPos());

    msc.parse("2I1D", 0);
    assertEquals(-1, msc.getStartPos());

    msc.parse("2I2D", 0);
    assertEquals(-1, msc.getStartPos());

    msc.parse("2I4D", 0);
    assertEquals(0, msc.getStartPos());

    msc.parse("42=6D1=2X53=", 0);
    assertEquals(3, msc.getStartPos());

//    179=1I3=1I11=1I4=7I1X8=3I3=2I5=2I22=1X6= picked ms/ns= 11/24739, would have been 18/24747
    msc.parse("179=1I3=1I11=1I4=7I1X8=3I3=2I5=2I22=1X6=", 24747);
    assertEquals(24739, msc.getStartPos());
  }

  public void testMaxShifts() {

    final MaxShiftCigarParser msc = new MaxShiftCigarParser();

    msc.parse("1I", 0);
    assertEquals(7, msc.getMaxShift());

    msc.parse("2I", 0);
    assertEquals(7, msc.getMaxShift());

    msc.parse("10I", 0);
    assertEquals(8, msc.getMaxShift()); //10 / 2 + 3

    msc.parse("9D", 0);
    assertEquals(7, msc.getMaxShift()); //10 / 2 + 3

    msc.parse("10D", 0);
    assertEquals(8, msc.getMaxShift()); //10 / 2 + 3

    msc.parse("2I1D", 0);
    assertEquals(7, msc.getMaxShift());

    msc.parse("6I6D", 0);
    assertEquals(7, msc.getMaxShift());

    msc.parse("10I12D", 0);
    assertEquals(9, msc.getMaxShift());

    msc.parse("10I22D", 0);
    assertEquals(14, msc.getMaxShift());

    msc.parse("42=6D1=2X53=", 0);
    assertEquals(7, msc.getMaxShift());

//    179=1I3=1I11=1I4=7I1X8=3I3=2I5=2I22=1X6= picked ms/ns= 11/24739, would have been 18/24747
    msc.parse("179=1I3=1I11=1I4=7I1X8=3I3=2I5=2I22=1X6=", 24747);
    assertEquals(11, msc.getMaxShift());
  }

  public void testSoftClip() {
    final MaxShiftCigarParser msc = new MaxShiftCigarParser();
    msc.parse("27S30=", 30);
    assertEquals(27, msc.getSoftClipStartOffset());
    assertEquals(30, msc.getStartPos());

    msc.parse("30=27S", 30);
    assertEquals(0, msc.getSoftClipStartOffset());
    assertEquals(30, msc.getStartPos());
  }
}
