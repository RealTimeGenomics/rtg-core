/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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

import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

import junit.framework.TestCase;

/**
 */
public class DiscordantReadSetTest extends TestCase {

  public void test() {
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0));
    assertEquals(42, drs.unionPosition());
    assertEquals(102, drs.flushPosition());
    final String exp = ""
      + "DiscordantReadSet:" + LS
      + "union=Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "intersection=Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "constraints count=1" + LS
      + "    Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      ;
    assertEquals(exp, drs.toString());
    //assertEquals("1 f 42 50 s 142 150 f 42 50 s 142 150 " + LS, drs.outputString().replace("\t", " "));
    assertEquals(1, drs.getCounts());
    assertEquals("Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0", drs.getIntersection().toString());
    assertEquals("Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0", drs.getUnion().toString());
    assertEquals("f", drs.getSequenceName());
  }

  public void testAdd() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);

    assertTrue(drs.belongs(bg));

    final BreakpointConstraint bg1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 40, 45, 140, 145, 180, 185), 0, 10.0);
    assertFalse(drs.belongs(bg1));

    drs.add(bg1);
    final String exp = ""
      + "DiscordantReadSet:" + LS
      + "union=Break-point constraint:UU x=40,50:f y=140,150:s r=180,198 Gap: mean=0.0 std.dev.=7.1" + LS
      + "intersection=null" + LS
      + "constraints count=2" + LS
      + "    Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "    Break-point constraint:UU x=40,45:f y=140,145:s r=180,185 Gap: mean=0.0 std.dev.=10.0" + LS
      ;
    assertEquals(exp, drs.toString());
    assertEquals(2, drs.getCounts());
    assertNull(drs.getIntersection());
    assertEquals("Break-point constraint:UU x=40,50:f y=140,150:s r=180,198 Gap: mean=0.0 std.dev.=7.1", drs.getUnion().toString());

    assertTrue(drs.belongs(bg));
    assertTrue(drs.belongs(bg1));

    final BreakpointGeometry bgTiny = new BreakpointGeometry(Orientation.UU, "f", "s", 42, 43, 143, 144, 185, 186);
    assertFalse(drs.belongs(bgTiny));
  }

  //There is an intersection but test case doesnt intersect it but does intersect one of the members
  public void testBelongs() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 188, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);

    assertTrue(drs.belongs(bg));

    final BreakpointConstraint bg1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 40, 45, 140, 149, 180, 189), 0, 10.0);
    assertTrue(drs.belongs(bg1));

    drs.add(bg1);
    final String exp = ""
      + "DiscordantReadSet:" + LS
      + "union=Break-point constraint:UU x=40,50:f y=140,150:s r=180,198 Gap: mean=0.0 std.dev.=7.1" + LS
      + "intersection=Break-point constraint:UU x=42,45:f y=143,147:s r=188,189 Gap: mean=0.0 std.dev.=7.1" + LS
      + "constraints count=2" + LS
      + "    Break-point constraint:UU x=42,50:f y=142,150:s r=188,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "    Break-point constraint:UU x=40,45:f y=140,149:s r=180,189 Gap: mean=0.0 std.dev.=10.0" + LS
      ;
    assertEquals(exp, drs.toString());

    assertTrue(drs.belongs(bg));
    assertTrue(drs.belongs(bg1));

    final BreakpointGeometry bgTiny1 = new BreakpointGeometry(Orientation.UU, "f", "s", 45, 47, 146, 148, 191, 193);
    assertTrue(drs.belongs(bgTiny1));

    final BreakpointGeometry bgTiny2 = new BreakpointGeometry(Orientation.UU, "f", "s", 43, 45, 144, 146, 187, 189);
    assertTrue(drs.belongs(bgTiny2));
  }

  public void testAddAll() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);

    final BreakpointConstraint bg1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 40, 45, 140, 145, 180, 185), 0, 10.0);
    final DiscordantReadSet drs1 = new DiscordantReadSet("f", 60, bg1);

    drs.addAll(drs1);
    final String exp = ""
      + "DiscordantReadSet:" + LS
      + "union=Break-point constraint:UU x=40,50:f y=140,150:s r=180,198 Gap: mean=0.0 std.dev.=7.1" + LS
      + "intersection=null" + LS
      + "constraints count=2" + LS
      + "    Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "    Break-point constraint:UU x=40,45:f y=140,145:s r=180,185 Gap: mean=0.0 std.dev.=10.0" + LS
      ;
    assertEquals(exp, drs.toString());

    final BreakpointConstraint bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 40, 45, 140, 145, 180, 185), 0, 10.0);
    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 60, bg2);
    drs2.addAll(drs);
    final String exp2 = ""
      + "DiscordantReadSet:" + LS
      + "union=Break-point constraint:UU x=40,50:f y=140,150:s r=180,198 Gap: mean=0.0 std.dev.=5.8" + LS
      + "intersection=null" + LS
      + "constraints count=3" + LS
      + "    Break-point constraint:UU x=40,45:f y=140,145:s r=180,185 Gap: mean=0.0 std.dev.=10.0" + LS
      + "    Break-point constraint:UU x=42,50:f y=142,150:s r=190,198 Gap: mean=0.0 std.dev.=10.0" + LS
      + "    Break-point constraint:UU x=40,45:f y=140,145:s r=180,185 Gap: mean=0.0 std.dev.=10.0" + LS
      ;
    assertEquals(exp2, drs2.toString());

  }

  public void testBelongsSameSequenceDirection() {
    final BreakpointGeometry bgFirst = new BreakpointGeometry(Orientation.UU, "f", "f", 45, 53, 50, 58, 95, 108);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, new BreakpointConstraint(bgFirst, 5, 10.0));
    final BreakpointGeometry bgSecond = new BreakpointGeometry(Orientation.UU, "f", "f", 50, 58, 45, 53, 95, 103);
    assertFalse(drs.belongs(bgSecond));
    final BreakpointGeometry bgThird = new BreakpointGeometry(Orientation.UU, "f", "f", 50, 58, 55, 60, 95, 97);
    assertTrue(drs.belongs(bgThird));


    // pair is at the same position
    final BreakpointGeometry bgSame = new BreakpointGeometry(Orientation.UU, "f", "f", 45, 53, 45, 53, 95, 100);
    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 60, new BreakpointConstraint(bgSame, 5, 10.0));
    assertTrue(drs.belongs(bgSame));
    assertFalse(drs2.belongs(bgSecond));

    // constraint is reversed
    final DiscordantReadSet drs3 = new DiscordantReadSet("f", 60, new BreakpointConstraint(bgSecond, 5, 10.0));
    assertFalse(drs3.belongs(bgFirst));

    final BreakpointGeometry bgFourth = new BreakpointGeometry(Orientation.UU, "f", "f", 51, 58, 45, 53, 95, 97);
    assertTrue(drs3.belongs(bgFourth));

  }
}
