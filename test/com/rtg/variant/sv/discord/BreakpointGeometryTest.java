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
import junit.framework.TestCase;

/**
 */
public class BreakpointGeometryTest extends TestCase {

  public void testConstructor() {
    final AbstractBreakpointGeometry bc = new BreakpointGeometry(Orientation.UU, "a", "b", 42, 151, 256, 365, 298, 407);
    bc.integrity();
    assertEquals(Orientation.UU, bc.getOrientation());
    assertEquals(42, bc.getXLo());
    assertEquals(256, bc.getYLo());

    assertEquals(/* 42 + 109 */151, bc.getXHi());
    assertEquals(/* 256 + 109 */ 365, bc.getYHi());

    assertEquals(298, bc.getRLo());
    assertEquals(407, bc.getRHi());

  }

  public void testConstructor2() {
    final AbstractBreakpointGeometry bc = new BreakpointGeometry(Orientation.UD, "a", "b", 42, 151, 256, 147, -214, -105);
    bc.integrity();
    assertEquals(Orientation.UD, bc.getOrientation());
    assertEquals(42, bc.getXLo());
    assertEquals(256, bc.getYLo());


    assertEquals(/* 42 + 109 */151, bc.getXHi());
    assertEquals(/* 256 - 109 */ 147, bc.getYHi());

    assertEquals(-214, bc.getRLo());
    assertEquals(-105, bc.getRHi());

  }

  public void testConstructor3() {
    final AbstractBreakpointGeometry bc = new BreakpointGeometry(Orientation.DD, "a", "b", 42, -67, 256, 147, -298, -189);
    bc.integrity();
    assertEquals(Orientation.DD, bc.getOrientation());
    assertEquals(42, bc.getXLo());
    assertEquals(256, bc.getYLo());

    assertEquals(/* 42 - 109 */ -67, bc.getXHi());
    assertEquals(/* 256 - 109 */ 147, bc.getYHi());

    assertEquals(-298, bc.getRLo());
    assertEquals(-189, bc.getRHi());
  }

  public void testConstructor4() {
    final AbstractBreakpointGeometry bc = new BreakpointGeometry(Orientation.DU, "a", "b", 256, 147, 42, 151, -214, -105);
    bc.integrity();
    assertEquals(Orientation.DU, bc.getOrientation());
    assertEquals("a", bc.getXName());
    assertEquals("b", bc.getYName());
    assertEquals(256, bc.getXLo());
    assertEquals(42, bc.getYLo());

    assertEquals(/* 256 - 109 */ 147, bc.getXHi());
    assertEquals(/* 42 + 109 */151, bc.getYHi());

    assertEquals(-214, bc.getRLo());
    assertEquals(-105, bc.getRHi());

    final String exp = "Break-point constraint:DU x=256,147:a y=42,151:b r=-214,-105";
    assertEquals(exp, bc.toString());

    final AbstractBreakpointGeometry fl = bc.flip();
    fl.integrity();
    assertEquals(Orientation.UD, fl.getOrientation());
    assertEquals("b", fl.getXName());
    assertEquals("a", fl.getYName());
    assertEquals(42, fl.getXLo());
    assertEquals(256, fl.getYLo());

    assertEquals(/* 42 + 109 */151, fl.getXHi());
    assertEquals(/* 256 - 109 */ 147, fl.getYHi());

    assertEquals(-214, fl.getRLo());
    assertEquals(-105, fl.getRHi());
  }

  public void testMakeGeometry1() {
    //see testConstructor4
    final BreakpointGeometry bg = BreakpointGeometry.makeGeometry("a", 256, 147, "b", 42, 151, -214, -105);
    final String exp = "Break-point constraint:DU x=256,147:a y=42,151:b r=-214,-105";
    assertEquals(exp, bg.toString());
  }

  public void testMakeGeometry2() {
    //see testConstructor4
    final BreakpointGeometry bg = BreakpointGeometry.makeGeometry("a", 147, 256, "b", 151, 42, 105, 214);
    final String exp = "Break-point constraint:UD x=147,256:a y=151,42:b r=105,214";
    assertEquals(exp, bg.toString());
  }

  //  public void testUnionBug() {
  //    final AbstractBreakpointGeometry bg = new BreakpointGeometry(Orientation.UU, "a", "a", 52875, 53251, 49976, 50352, 102877, 103227);
  //    bg.globalIntegrity();
  //    //System.err.println("bg");
  //    System.err.println(bg.gnuPlot());
  //
  //    final AbstractBreakpointGeometry bu = new BreakpointGeometry(Orientation.UU, "a", "a", 49913, 50362, 52815, 53232, 102774, 103218).flip();
  //    bu.globalIntegrity();
  //    //System.err.println("bu");
  //    System.err.println(bu.gnuPlot());
  //
  //    //System.err.println("in");
  //    System.err.println(bg.intersect(bu).gnuPlot());
  //
  //    //System.err.println("un");
  //    System.err.println(bg.union(bu).gnuPlot());
  //  }
}
