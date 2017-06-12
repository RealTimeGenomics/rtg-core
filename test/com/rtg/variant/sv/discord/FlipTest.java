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
    protected static final int START_FIRST = 42;
    protected static final int LENGTH = 7;
    protected static final int START_SECOND = 256;
    protected static final int START_END_EX = START_FIRST + LENGTH;
    protected static final ReadGroupStats RGS = new MockReadGroupStats(100);
    protected static final int MAX_GAP = BreakpointConstraint.gapMax(RGS);
    protected static final int MIN_GAP = BreakpointConstraint.gapMin(RGS);

    protected BreakpointConstraint mBcFirst;

    TestConstructor(final int firstFlag) {
      assertEquals(113, MAX_GAP);
      assertEquals(87, MIN_GAP);
      final MockSam recFirst = new MockSam();
      recFirst.setAlignmentStart(START_FIRST);
      recFirst.setAlignmentEnd(START_END_EX - 1);
      recFirst.setMateAlignmentStart(START_SECOND);
      recFirst.setFlags(firstFlag);
      mBcFirst = new BreakpointConstraint(recFirst, MO, RGS);
      mBcFirst.integrity();

      checkFirst();

      final MockSam recSecond = new MockSam();
      recSecond.setAlignmentStart(START_SECOND);
      recSecond.setAlignmentEnd(START_SECOND + LENGTH - 1);
      recSecond.setMateAlignmentStart(START_FIRST);
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
        assertEquals(START_END_EX, mBcFirst.getX());
        assertEquals(START_SECOND + LENGTH, mBcFirst.getY());
        assertEquals(START_END_EX + MAX_GAP, mBcFirst.getZ());
        assertEquals(START_SECOND + LENGTH + MAX_GAP, mBcFirst.getW());
        assertEquals(START_END_EX + START_SECOND + LENGTH + MIN_GAP, mBcFirst.getR());
        assertEquals(START_END_EX + START_SECOND + LENGTH + MAX_GAP, mBcFirst.getS());
      }
    };
  }

  public void testConstructor2() {
    new TestConstructor(OrientationTest.FLAG_FIRST_FR) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.UD, mBcFirst.getOrientation());
        assertEquals(START_END_EX, mBcFirst.getX());
        assertEquals(START_SECOND - 1, mBcFirst.getY());
        assertEquals(START_END_EX + MAX_GAP, mBcFirst.getZ());
        assertEquals(START_SECOND - 1 - MAX_GAP, mBcFirst.getW());
        assertEquals(START_END_EX - (START_SECOND - 1) + MIN_GAP, mBcFirst.getR());
        assertEquals(START_END_EX - (START_SECOND - 1) + MAX_GAP, mBcFirst.getS());
      }
    };
  }

  public void testConstructor3() {
    new TestConstructor(OrientationTest.FLAG_FIRST_RR) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.DD, mBcFirst.getOrientation());
        assertEquals(START_FIRST - 1, mBcFirst.getX());
        assertEquals(START_SECOND - 1, mBcFirst.getY());
        assertEquals(START_FIRST - 1 - MAX_GAP, mBcFirst.getZ());
        assertEquals(START_SECOND - 1 - MAX_GAP, mBcFirst.getW());
        assertEquals(-START_FIRST + 1 - START_SECOND + 1 + MIN_GAP, mBcFirst.getR());
        assertEquals(-START_FIRST + 1 - START_SECOND + 1 + MAX_GAP, mBcFirst.getS());
      }
    };
  }

  public void testConstructor4() {
    new TestConstructor(OrientationTest.FLAG_FIRST_RF) {

      @Override
      void checkFirst() {
        assertEquals(Orientation.DU, mBcFirst.getOrientation());
        assertEquals(START_FIRST - 1, mBcFirst.getX());
        assertEquals(START_SECOND + LENGTH, mBcFirst.getY());
        assertEquals(START_FIRST - 1 - MAX_GAP, mBcFirst.getZ());
        assertEquals(START_SECOND + LENGTH + MAX_GAP, mBcFirst.getW());
        assertEquals(-START_FIRST + 1 + START_SECOND + LENGTH + MIN_GAP, mBcFirst.getR());
        assertEquals(-START_FIRST + 1 + START_SECOND + LENGTH + MAX_GAP, mBcFirst.getS());
      }
    };
  }

  public void testFlip() {
    new TestConstructor(OrientationTest.FLAG_FIRST_FR) {

      @Override
      void checkFirst() {
        assertEquals(mBcFirst.getX(), mBcFirst.flip().getY());
        assertEquals(mBcFirst.getY(), mBcFirst.flip().getX());
        assertEquals(mBcFirst.getZ(), mBcFirst.flip().getW());
        assertEquals(mBcFirst.getW(), mBcFirst.flip().getZ());
        assertEquals(mBcFirst.getR(), mBcFirst.flip().getR());
        assertEquals(mBcFirst.getS(), mBcFirst.flip().getS());
        assertEquals(Orientation.DU, mBcFirst.flip().getOrientation());
        assertEquals(mBcFirst, mBcFirst.flip().flip());
      }
    };
  }
}
