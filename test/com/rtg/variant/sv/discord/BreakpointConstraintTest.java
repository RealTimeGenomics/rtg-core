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

import static com.rtg.util.StringUtils.LS;

import com.rtg.sam.SamUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.variant.sv.ReadGroupStats;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class BreakpointConstraintTest extends TestCase {

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(BreakpointConstraintTest.class);
    suite.addTestSuite(FlipTest.class);
    return suite;
  }

  private static class MockReadGroupStats extends ReadGroupStats {
    static final double STD_DEV =  3.14159;
    private final double mMean;
    MockReadGroupStats(double mean) {
      super("id", 1000);
      mMean = mean;
    }
    @Override
    public double gapStdDev() {
      return STD_DEV;
    }
    @Override
    public double gapMean() {
      return mMean;
    }
  }

  private static BreakpointConstraint makeConstraint(BreakpointGeometry bg) {
    return new BreakpointConstraint(bg, 100, MockReadGroupStats.STD_DEV);
  }
  private static BreakpointConstraint makeConstraint(String xSeq, int x, boolean xUp, String ySeq, int y, boolean yUp, double mean) {
    final ReadGroupStats rgs = new MockReadGroupStats(mean);
    final MockSam rec = new MockSam();
    rec.setReferenceName(xSeq);
    rec.setAlignmentStart(x);
    rec.setAlignmentEnd(x + 7);
    rec.setReadPairedFlag(true);
    rec.setMateReferenceName(ySeq);
    rec.setMateAlignmentStart(y);
    rec.setAttribute(SamUtils.ATTRIBUTE_MATE_END, y + 7);
    rec.setFirstOfPairFlag(true);
    rec.setReadNegativeStrandFlag(!xUp);
    rec.setMateNegativeStrandFlag(!yUp);

    final MachineOrientation mo = MachineOrientation.FR;
    final BreakpointConstraint bc = new BreakpointConstraint(rec, mo, rgs);
    bc.integrity();
    return bc;
  }

  public void test() {
    final BreakpointConstraint bc = makeConstraint("x", 110, true, "y", 110, true, 10);
    final String exp = ""
        + "Break-point constraint:UU x=118,141:x y=118,141:y r=236,259 Gap: mean=246.0 std.dev.=3.1"
        ;
    assertEquals(exp, bc.toString());
    assertEquals(246.0, bc.rMean());
    assertEquals(3.14159, bc.rStdDev());
    assertEquals("BreakpointPosition [lo=118, position=123, hi=128, positionY=123]", bc.position().toString());
    assertEquals("BreakpointPosition [lo=118, position=123, hi=128, positionY=123]", bc.flip().position().toString());
    final String gpx = ""
        + "118.0 118.0" + LS
        + "118.0 118.0" + LS
        + "118.0 141.0" + LS
        + "118.0 141.0" + LS
        + "141.0 118.0" + LS
        + "141.0 118.0" + LS
        + "118.0 118.0" + LS
        ;
    assertEquals(gpx, bc.gnuPlot().replace("\t ", " "));
    final String gpmx = ""
        + "118.0 128.0" + LS
        + "128.0 118.0" + LS
        ;
    assertEquals(gpmx, bc.gnuPlotMean().replace("\t ", " "));

    assertNull(bc.outputGeometry("z"));
    assertTrue(bc == bc.outputGeometry("x"));
    assertEquals("y", bc.outputGeometry("y").getXName());
  }

  public void testPositionDD1() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -55.6, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(15, p.lo());
    assertEquals(17, p.position());
    assertEquals(20, p.hi());
    assertEquals(37, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(35, q.lo());
    assertEquals(37, q.position());
    assertEquals(40, q.hi());
    assertEquals(17, q.positionAlt());
  }

  public void testPositionDD2() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -55.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(15, p.lo());
    assertEquals(17, p.position());
    assertEquals(20, p.hi());
    assertEquals(37, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(35, q.lo());
    assertEquals(37, q.position());
    assertEquals(40, q.hi());
    assertEquals(17, q.positionAlt());
  }

  public void testPositionDD3() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -53.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(13, p.lo());
    assertEquals(16, p.position());
    assertEquals(20, p.hi());
    assertEquals(36, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(33, q.lo());
    assertEquals(36, q.position());
    assertEquals(40, q.hi());
    assertEquals(16, q.positionAlt());
  }

  public void testPositionDD4() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -50.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionDD5() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -47.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(13, p.position());
    assertEquals(17, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(37, q.hi());
    assertEquals(13, q.positionAlt());
  }

  public void testPositionDD6() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -45.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(13, p.position());
    assertEquals(16, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(36, q.hi());
    assertEquals(13, q.positionAlt());
  }

  public void testPositionDD7() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "x", "y", 20, 10, 40, 30, -55, -45), -44.4, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(13, p.position());
    assertEquals(16, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(36, q.hi());
    assertEquals(13, q.positionAlt());
  }

  public void testPositionDU1() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 12.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(14, p.lo());
    assertEquals(17, p.position());
    assertEquals(20, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(36, q.hi());
    assertEquals(17, q.positionAlt());
  }

  public void testPositionDU2() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 16.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(14, p.lo());
    assertEquals(17, p.position());
    assertEquals(20, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(36, q.hi());
    assertEquals(17, q.positionAlt());
  }

  public void testPositionDU3() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 16.6, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(13, p.lo());
    assertEquals(16, p.position());
    assertEquals(20, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(37, q.hi());
    assertEquals(16, q.positionAlt());
  }

  public void testPositionDU4() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 20.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionDU5() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 25.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(12, p.position());
    assertEquals(15, p.hi());
    assertEquals(37, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(35, q.lo());
    assertEquals(37, q.position());
    assertEquals(40, q.hi());
    assertEquals(12, q.positionAlt());
  }

  public void testPositionDU6() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 26.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(12, p.position());
    assertEquals(14, p.hi());
    assertEquals(38, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(36, q.lo());
    assertEquals(38, q.position());
    assertEquals(40, q.hi());
    assertEquals(12, q.positionAlt());
  }

  public void testPositionDU7() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.DU, "x", "y", 20, 10, 30, 40, 16, 27), 27.6, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(12, p.position());
    assertEquals(14, p.hi());
    assertEquals(38, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(36, q.lo());
    assertEquals(38, q.position());
    assertEquals(40, q.hi());
    assertEquals(12, q.positionAlt());
  }

  public void testPositionUD1() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -25.5, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(14, p.position());
    assertEquals(18, p.hi());
    assertEquals(36, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(32, q.lo());
    assertEquals(36, q.position());
    assertEquals(40, q.hi());
    assertEquals(14, q.positionAlt());
  }

  public void testPositionUD2() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -22.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(14, p.position());
    assertEquals(18, p.hi());
    assertEquals(36, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(32, q.lo());
    assertEquals(36, q.position());
    assertEquals(40, q.hi());
    assertEquals(14, q.positionAlt());
  }

  public void testPositionUD3() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -21.4, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(14, p.position());
    assertEquals(19, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(31, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(14, q.positionAlt());
  }

  public void testPositionUD4() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -20.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionUD5() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -15.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(15, p.lo());
    assertEquals(17, p.position());
    assertEquals(20, p.hi());
    assertEquals(32, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(32, q.position());
    assertEquals(35, q.hi());
    assertEquals(17, q.positionAlt());
  }

  public void testPositionUD6() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -12.4, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(17, p.lo());
    assertEquals(18, p.position());
    assertEquals(20, p.hi());
    assertEquals(31, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(31, q.position());
    assertEquals(33, q.hi());
    assertEquals(18, q.positionAlt());
  }

  public void testPositionUD7() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 10, 20, 40, 30, -22, -12), -12.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(17, p.lo());
    assertEquals(18, p.position());
    assertEquals(20, p.hi());
    assertEquals(31, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(31, q.position());
    assertEquals(33, q.hi());
    assertEquals(18, q.positionAlt());
  }

  public void testPositionUU1() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 41.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(11, p.position());
    assertEquals(12, p.hi());
    assertEquals(31, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(31, q.position());
    assertEquals(32, q.hi());
    assertEquals(11, q.positionAlt());
  }

  public void testPositionUU2() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 42.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(11, p.position());
    assertEquals(12, p.hi());
    assertEquals(31, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(31, q.position());
    assertEquals(32, q.hi());
    assertEquals(11, q.positionAlt());
  }

  public void testPositionUU3() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 45.5, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(13, p.position());
    assertEquals(16, p.hi());
    assertEquals(33, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(33, q.position());
    assertEquals(36, q.hi());
    assertEquals(13, q.positionAlt());
  }

  public void testPositionUU4() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 50.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(10, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(30, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionUU5() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 51.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(11, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(31, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionUU6() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 52.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(11, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(31, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testPositionUU7() {
    final BreakpointConstraint bc = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 10, 20, 30, 40, 42, 52), 53.0, 12.3);
    //displayPosition(bc);
    final BreakpointPosition p = bc.position();
    assertEquals(11, p.lo());
    assertEquals(15, p.position());
    assertEquals(20, p.hi());
    assertEquals(35, p.positionAlt());

    final BreakpointPosition q = bc.flip().position();
    assertEquals(31, q.lo());
    assertEquals(35, q.position());
    assertEquals(40, q.hi());
    assertEquals(15, q.positionAlt());
  }

  public void testMakeConstraint() {
    assertEquals(Orientation.UU, makeConstraint("x", 110, true, "y", 110, true, 10).getOrientation());
    assertEquals(Orientation.UD, makeConstraint("x", 110, true, "y", 110, false, 10).getOrientation());
    assertEquals(Orientation.DU, makeConstraint("x", 110, false, "y", 110, true, 10).getOrientation());
    assertEquals(Orientation.DD, makeConstraint("x", 110, false, "y", 110, false, 10).getOrientation());
  }

  public void testOverlapUU() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint x = makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 141, 118, 141, 236, 259));
    checkOverlap(bc1, bc2, x);
  }

  public void testOverlapUD() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, false, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, true, "y", 110, false, 10);
    final BreakpointConstraint x = makeConstraint(new BreakpointGeometry(Orientation.UD, "x", "y", 118, 141, 109, 86, 9, 32));
    checkOverlap(bc1, bc2, x);
  }

  public void testOverlapSeq1() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, true, "z", 110, true, 10);
    checkOverlap(bc1, bc2, null);
  }

  public void testOverlapSeq2() {
    final BreakpointConstraint bc1 = makeConstraint("z", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, true, "y", 110, true, 10);
    checkOverlap(bc1, bc2, null);
  }

  public void testOverlapOrientation1() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, false, "y", 110, true, 10);
    checkOverlap(bc1, bc2, null);
  }

  public void testOverlapOrientation2() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, true, "y", 110, false, 10);
    checkOverlap(bc1, bc2, null);
  }

  public void testOverlapOrientation3() {
    final BreakpointConstraint bc1 = makeConstraint("x", 110, true, "y", 110, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 110, false, "y", 110, false, 10);
    checkOverlap(bc1, bc2, null);
  }

  //see BreakpointConstraintTest.xslx
  public void testOverlapStatistics() {
    final AbstractBreakpointGeometry exp = new BreakpointGeometry(Orientation.UU, "x", "y", 118, 139, 118, 139, 236, 257);
    final BreakpointConstraint bc1 = makeConstraint("x", 100 + 10, true, "y", 100 + 10, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 100 + 9, true, "y", 100 + 9, true, 10);
    //System.err.println("bc1=" + bc1);
    //System.err.println("bc2=" + bc2);
    bc1.integrity();
    bc2.integrity();
    final BreakpointConstraint i1 = bc1.intersect(bc2);
    //System.err.println(" i1=" + i1);
    assertEquals(exp, i1);
    checkStats(i1);
    final BreakpointConstraint i2 = bc2.intersect(bc1);
    //System.err.println(" i2=" + i2);
    assertEquals(exp, i2);
    checkStats(i2);
  }

  private void checkStats(BreakpointConstraint bc) {
    assertEquals(245.0, bc.rMean(), 0.000001);
    assertEquals(2.221440, bc.rStdDev(), 0.000001);
  }

  private static void checkOverlap(boolean xUp, boolean yUp, int x0, int y0, int x1, int y1, BreakpointConstraint exp) {
    final BreakpointConstraint bc1 = makeConstraint("x", 100 + x0, xUp, "y", 100 + y0, yUp, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 100 + x1, xUp, "y", 100 + y1, yUp, 10);
    checkOverlap(bc1, bc2, exp);
  }

  public void testOverlapDiagonal() {
    checkOverlap(true, true, 10, 10, 19, 24, null);
    checkOverlap(true, true, 10, 10, 19, 23, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 127, 128, 131, 132, 258, 259)));
  }

  public void testOverlapXYLimit() {
    checkOverlap(true, true, 10, 10, 33, 10, null);
    checkOverlap(true, true, 10, 10, 32, 10, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 140, 141, 118, 119, 258, 259)));
  }

  public void testOverlapCorner() {
    checkOverlap(true, true, 10, 10, -8, 33, null);
    checkOverlap(true, true, 10, 10, -8, 32, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 119, 140, 141, 258, 259)));
  }

  public void testOverlapContaining() {
    checkOverlap(true, true, 10, 10, 9, 9, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 139, 118, 139, 236, 257)));
  }

  public void testOverlapMidpointXY() {
    checkOverlap(true, true, 10, 10, 15, -13, null);
    checkOverlap(true, true, 10, 10, 15, -12, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 123, 124, 118, 119, 241, 242)));
  }

  /**
   * @param bc1 first breakpoint
   * @param bc2 second breakpoint
   */
  private static <T extends AbstractBreakpointGeometry> void checkOverlap(T bc1, T bc2, T expected) {
    //System.err.println("bc1=" + bc1);
    //System.err.println("bc2=" + bc2);
    //System.err.println("exp=" + expected);
    bc1.integrity();
    bc2.integrity();
    assertTrue(bc1.overlap(bc1));
    assertTrue(bc2.overlap(bc2));
    assertEquals(expected != null, bc1.overlap(bc2));
    assertEquals(expected != null, bc2.overlap(bc1));
    chOver(expected, bc1, bc2);
    chOver(expected, bc1.flip(), bc2.flip());
    chOver(expected, bc2, bc1);
    chOver(expected, bc2.flip(), bc1.flip());
  }

  private static <T extends AbstractBreakpointGeometry> void chOver(T expected, T bc1, T bc2) {
    final AbstractBreakpointGeometry i1 = bc1.intersect(bc2);
    assertTrue(i1 == null || i1 instanceof BreakpointConstraint);
    assertEquals(expected, i1);
  }

  /**
   * @param bc1 first breakpoint
   * @param bc2 second breakpoint
   */
  private static <T extends AbstractBreakpointGeometry> void checkUnion(T bc1, T bc2, T expected) {
    bc1.integrity();
    bc2.integrity();
    assertEquals(bc1, bc1.union(bc1));
    assertEquals(bc2, bc2.union(bc2));
    //System.err.println("bc1=" + bc1);
    //System.err.println("bc2=" + bc2);
    //System.err.println("exp=" + expected);
    final AbstractBreakpointGeometry i1 = bc1.union(bc2);
    assertEquals(expected, i1);
    assertEquals(expected, bc2.union(bc1));
  }

  private static void checkUnion(boolean xUp, boolean yUp, int x0, int y0, int x1, int y1, BreakpointConstraint exp) {
    checkUnion("x", "y", xUp, yUp, x0, y0, x1, y1, exp);
  }

  private static void checkUnion(String x, String y, boolean xUp, boolean yUp, int x0, int y0, int x1, int y1, BreakpointConstraint exp) {
    final BreakpointConstraint bc1 = makeConstraint(x, 100 + x0, xUp, y, 100 + y0, yUp, 10);
    final BreakpointConstraint bc2 = makeConstraint(x, 100 + x1, xUp, y, 100 + y1, yUp, 10);
    checkUnion(bc1, bc2, exp);
  }

  public void testUnionNull0() {
    final BreakpointConstraint bc1 = makeConstraint("a", 100, true, "b", 100, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("x", 100, true, "b", 100, true, 10);
    checkUnion(bc1, bc2, null);
  }

  public void testUnionNull1() {
    final BreakpointConstraint bc1 = makeConstraint("a", 100, true, "b", 100, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("a", 100, true, "x", 100, true, 10);
    checkUnion(bc1, bc2, null);
  }

  public void testUnionNull2() {
    final BreakpointConstraint bc1 = makeConstraint("a", 100, true, "b", 100, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("a", 100, true, "b", 100, false, 10);
    checkUnion(bc1, bc2, null);
  }

  public void testUnionNull3() {
    final BreakpointConstraint bc1 = makeConstraint("a", 100, true, "b", 100, true, 10);
    final BreakpointConstraint bc2 = makeConstraint("a", 100, false, "b", 100, true, 10);
    checkUnion(bc1, bc2, null);
  }

  public void testUnionDiagonal() {
    checkUnion(true, true, 10, 10, 19, 19, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 150, 118, 150, 236, 277)));
  }

  public void testUnionXYLimit() {
    checkUnion(true, true, 10, 10, 28, 10, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 159, 118, 141, 236, 277)));
  }

  public void testUnionCorner() {
    checkUnion(true, true, 10, 10, -8, 28, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 100, 141, 118, 159, 236, 259)));
  }

  public void testUnionContaining() {
    checkUnion(true, true, 10, 10, 9, 9, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 117, 141, 117, 141, 234, 259)));
  }

  public void testUnionMidpointXY() {
    checkUnion(true, true, 10, 10, 15, -8, makeConstraint(new BreakpointGeometry(Orientation.UU, "x", "y", 118, 146, 100, 141, 223, 259)));
  }

  public void testIsConcordant() {
    assertFalse(makeConstraint("x", 110, true, "y", 110 + 9, false, 10).isConcordant());
    assertFalse(makeConstraint("x", 110, false, "x", 110 + 9, false, 10).isConcordant());
    assertFalse(makeConstraint("x", 110, true, "x", 110 + 9, true, 10).isConcordant());
  }

  private void checkIsConcordantUD(boolean expected, int offset) {
    final int mean = 100;
    final int basex = 42;
    final int length = 7;
    final int basey = basex + length + mean + 1;
    assertEquals(true, makeConstraint("x", basex, true, "x", basey, false, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basex, true, "x", basey + offset, false, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basex, true, "x", basey - offset, false, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basex + offset, true, "x", basey, false, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basex - offset, true, "x", basey, false, mean).isConcordant());
  }

  public void testIsConcordantUD() {
    checkIsConcordantUD(true, 13);
    checkIsConcordantUD(false, 14);
  }

  private void checkIsConcordantDU(boolean expected, int offset) {
    final int mean = 100;
    final int basex = 42;
    final int length = 7;
    final int basey = basex + length + mean + 1;
    assertEquals(true, makeConstraint("x", basey, false, "x", basex, true, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basey, false, "x", basex + offset, true, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basey, false, "x", basex - offset, true, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basey + offset, false, "x", basex, true, mean).isConcordant());
    assertEquals(expected, makeConstraint("x", basey - offset, false, "x", basex, true, mean).isConcordant());
  }

  public void testIsConcordantDU() {
    checkIsConcordantDU(true, 13);
    checkIsConcordantDU(false, 14);
  }

  private static BreakpointConstraint makeMateEnd(int y, final Integer mateEnd) {
    final ReadGroupStats rgs = new MockReadGroupStats(10.0);
    final MockSam rec = new MockSam();
    rec.setReferenceName("x");
    rec.setAlignmentStart(10);
    rec.setAlignmentEnd(10 + 7);
    rec.setReadPairedFlag(true);
    rec.setMateReferenceName("y");
    rec.setMateAlignmentStart(y);
    if (mateEnd != null) {
      rec.setAttribute(SamUtils.ATTRIBUTE_MATE_END, mateEnd);
    }
    rec.setFirstOfPairFlag(true);
    rec.setReadNegativeStrandFlag(true);
    rec.setMateNegativeStrandFlag(false);

    final MachineOrientation mo = MachineOrientation.FR;
    final BreakpointConstraint bc = new BreakpointConstraint(rec, mo, rgs);
    bc.integrity();
    return bc;
  }

  public void testMateEnd() {
    assertEquals(28, makeMateEnd(20, 27).getY());
    assertEquals(28, makeMateEnd(20, null).getY());
    assertEquals(31, makeMateEnd(20, 30).getY());
  }

  private static BreakpointConstraint makeAS(final boolean first, final Integer as, final Integer mateAS, boolean xUp, boolean yUp) {
    final ReadGroupStats rgs = new MockReadGroupStats(10.0);
    final MockSam rec = new MockSam();
    rec.setReferenceName("x");
    rec.setAlignmentStart(10);
    rec.setAlignmentEnd(10 + 7);
    rec.setReadPairedFlag(true);
    rec.setMateReferenceName("y");
    rec.setMateAlignmentStart(20);
    rec.setAttribute(SamUtils.ATTRIBUTE_MATE_END, 30);
    if (as != null) {
      rec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, as);
    }
    if (mateAS != null) {
      rec.setAttribute(SamUtils.ATTRIBUTE_MATE_ALIGNMENT_SCORE, mateAS);
    }
    rec.setFirstOfPairFlag(first);
    rec.setReadNegativeStrandFlag(!xUp);
    rec.setMateNegativeStrandFlag(!yUp);
    final MachineOrientation mo = MachineOrientation.FR;
    final BreakpointConstraint bc = new BreakpointConstraint(rec, mo, rgs);
    bc.integrity();
    return bc;
  }

  public void testAs() {
    final BreakpointConstraint bc = makeAS(true, null, null, false, true);
    assertEquals(9, bc.getX());
    assertEquals(31, bc.getY());
    final BreakpointConstraint bf = makeAS(false, null, null, false, true);
    assertEquals(9, bf.getX());
    assertEquals(31, bf.getY());
  }

  public void testAs0() {
    final BreakpointConstraint bc = makeAS(true, 0, 0, false, true);
    assertEquals(9, bc.getX());
    assertEquals(31, bc.getY());
    final BreakpointConstraint bf = makeAS(false, 0, 0, false, true);
    assertEquals(9, bf.getX());
    assertEquals(31, bf.getY());
  }

  public void testAs13UU() {
    final BreakpointConstraint bc = makeAS(true, 1, 3, true, true);
    assertEquals(17, bc.getX());
    assertEquals(28, bc.getY());
    final BreakpointConstraint bf = makeAS(false, 1, 3, true, true);
    assertEquals(28, bf.getY());
    assertEquals(17, bf.getX());
  }

  public void testAs13UD() {
    final BreakpointConstraint bc = makeAS(true, 1, 3, true, false);
    assertEquals(17, bc.getX());
    assertEquals(22, bc.getY());
    final BreakpointConstraint bf = makeAS(false, 1, 3, true, false);
    assertEquals(22, bf.getY());
    assertEquals(17, bf.getX());
  }

  public void testAs13DU() {
    final BreakpointConstraint bc = makeAS(true, 1, 3, false, true);
    assertEquals(10, bc.getX());
    assertEquals(28, bc.getY());
    final BreakpointConstraint bf = makeAS(false, 1, 3, false, true);
    assertEquals(28, bf.getY());
    assertEquals(10, bf.getX());
  }

  public void testAs13DD() {
    final BreakpointConstraint bc = makeAS(true, 1, 3, false, false);
    assertEquals(10, bc.getX());
    assertEquals(22, bc.getY());
    final BreakpointConstraint bf = makeAS(false, 1, 3, false, false);
    assertEquals(22, bf.getY());
    assertEquals(10, bf.getX());
  }
}
