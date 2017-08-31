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

package com.rtg.variant.sv.discord;

import com.rtg.util.machine.MachineOrientation;
import com.rtg.variant.sv.ReadGroupStats;

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
