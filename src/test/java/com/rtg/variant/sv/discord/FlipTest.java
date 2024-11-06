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

package com.rtg.variant.sv.discord;

import com.rtg.util.machine.MachineOrientation;
import com.rtg.variant.sv.ReadGroupStats;
import com.rtg.variant.sv.bndeval.Orientation;
import com.rtg.variant.sv.bndeval.OrientationTest;

import junit.framework.TestCase;

/**
 * Tests both BreakpointConstraint and FlippedProxyBreakpointConstraint.
 */
public class FlipTest extends TestCase {

  private static class MockReadGroupStats extends ReadGroupStats {
    private final double mMean;
    MockReadGroupStats(double mean) {
      super("id", 1000);
      mMean = mean;
    }
    @Override
    public double fragmentStdDev() {
      return Math.PI;
    }
    @Override
    public double gapMean() {
      return mMean;
    }
  }

  private abstract static class TestConstructor {

    protected static final MachineOrientation MO = MachineOrientation.FR;
    protected static final int LENGTH = 7;
    protected static final int START_FIRST = 42; // These constants are 0-based
    protected static final int START_SECOND = 256;
    protected static final int END_FIRST_EXCL = START_FIRST + LENGTH;
    protected static final int END_SECOND_EXCL = START_SECOND + LENGTH;
    protected static final ReadGroupStats RGS = new MockReadGroupStats(100);
    static {
      RGS.setNumDeviations(4.0);
    }
    protected static final int MAX_GAP = RGS.gapMax();
    protected static final int MIN_GAP = RGS.gapMin();

    protected BreakpointConstraint mBcFirst;

    TestConstructor(final int firstFlag) {
      assertEquals(113, MAX_GAP);
      assertEquals(87, MIN_GAP);
      final MockSam recFirst = new MockSam();
      recFirst.setAlignmentStart(START_FIRST + 1); // Various +1's here and below since htsjdk is 1-based
      recFirst.setAlignmentEnd(END_FIRST_EXCL - 1 + 1);  // And their alignmentEnd is inclusive
      recFirst.setMateAlignmentStart(START_SECOND + 1);
      recFirst.setFlags(firstFlag);
      mBcFirst = new BreakpointConstraint(recFirst, MO, RGS);
      mBcFirst.integrity();

      checkFirst();

      final MockSam recSecond = new MockSam();
      recSecond.setAlignmentStart(START_SECOND + 1);
      recSecond.setAlignmentEnd(END_SECOND_EXCL - 1 + 1);
      recSecond.setMateAlignmentStart(START_FIRST + 1);
      recSecond.setFlags(OrientationTest.flip(firstFlag));
      final BreakpointConstraint bcSecond = new BreakpointConstraint(recSecond, MO, RGS);
      assertEquals(mBcFirst, bcSecond);
    }
    abstract void checkFirst();
  }

  public void testConstructor1() {
    new TestConstructor(OrientationTest.FLAG_FIRST_FF) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.UU, mBcFirst.getOrientation());
        assertEquals(END_FIRST_EXCL, mBcFirst.getXLo());
        assertEquals(END_SECOND_EXCL, mBcFirst.getYLo());
        assertEquals(END_FIRST_EXCL + MAX_GAP, mBcFirst.getXHi());
        assertEquals(END_SECOND_EXCL + MAX_GAP, mBcFirst.getYHi());
        assertEquals(END_FIRST_EXCL + END_SECOND_EXCL + MIN_GAP, mBcFirst.getRLo());
        assertEquals(END_FIRST_EXCL + END_SECOND_EXCL + MAX_GAP, mBcFirst.getRHi());
      }
    };
  }

  public void testConstructor2() {
    new TestConstructor(OrientationTest.FLAG_FIRST_FR) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.UD, mBcFirst.getOrientation());
        assertEquals(END_FIRST_EXCL, mBcFirst.getXLo());
        assertEquals(START_SECOND, mBcFirst.getYLo());
        assertEquals(END_FIRST_EXCL + MAX_GAP, mBcFirst.getXHi());
        assertEquals(START_SECOND - MAX_GAP, mBcFirst.getYHi());
        assertEquals(END_FIRST_EXCL - START_SECOND + MIN_GAP, mBcFirst.getRLo());
        assertEquals(END_FIRST_EXCL - START_SECOND + MAX_GAP, mBcFirst.getRHi());
      }
    };
  }

  public void testConstructor3() {
    new TestConstructor(OrientationTest.FLAG_FIRST_RR) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.DD, mBcFirst.getOrientation());
        assertEquals(START_FIRST, mBcFirst.getXLo());
        assertEquals(START_SECOND, mBcFirst.getYLo());
        assertEquals(START_FIRST - MAX_GAP, mBcFirst.getXHi());
        assertEquals(START_SECOND - MAX_GAP, mBcFirst.getYHi());
        assertEquals(-START_FIRST - START_SECOND + MIN_GAP, mBcFirst.getRLo());
        assertEquals(-START_FIRST - START_SECOND + MAX_GAP, mBcFirst.getRHi());
      }
    };
  }

  public void testConstructor4() {
    new TestConstructor(OrientationTest.FLAG_FIRST_RF) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.DU, mBcFirst.getOrientation());
        assertEquals(START_FIRST, mBcFirst.getXLo());
        assertEquals(END_SECOND_EXCL, mBcFirst.getYLo());
        assertEquals(START_FIRST - MAX_GAP, mBcFirst.getXHi());
        assertEquals(END_SECOND_EXCL + MAX_GAP, mBcFirst.getYHi());
        assertEquals(-START_FIRST + END_SECOND_EXCL + MIN_GAP, mBcFirst.getRLo());
        assertEquals(-START_FIRST + END_SECOND_EXCL + MAX_GAP, mBcFirst.getRHi());
      }
    };
  }

  public void testFlip() {
    new TestConstructor(OrientationTest.FLAG_FIRST_FR) {

      @Override
      void checkFirst() {
        assertEquals(mBcFirst.getXLo(), mBcFirst.flip().getYLo());
        assertEquals(mBcFirst.getYLo(), mBcFirst.flip().getXLo());
        assertEquals(mBcFirst.getXHi(), mBcFirst.flip().getYHi());
        assertEquals(mBcFirst.getYHi(), mBcFirst.flip().getXHi());
        assertEquals(mBcFirst.getRLo(), mBcFirst.flip().getRLo());
        assertEquals(mBcFirst.getRHi(), mBcFirst.flip().getRHi());
        assertEquals(Orientation.DU, mBcFirst.flip().getOrientation());
        assertEquals(mBcFirst, mBcFirst.flip().flip());
      }
    };
  }
}
