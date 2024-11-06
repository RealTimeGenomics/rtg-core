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
