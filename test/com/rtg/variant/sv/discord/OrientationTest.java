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

import com.rtg.util.TestUtils;
import com.rtg.util.machine.MachineOrientation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class OrientationTest extends TestCase {
  // FLAG_{first/second in pair}_{read mapped to forward or reverse}{mate mapped to forward or reverse}
  static final int FLAG_FIRST_FR = 97;
  static final int FLAG_SECOND_RF = 145;

  static final int FLAG_FIRST_RF = 81;
  static final int FLAG_SECOND_FR = 161;

  static final int FLAG_FIRST_RR = 113;
  static final int FLAG_SECOND_RR = 177;

  static final int FLAG_FIRST_FF = 65;
  static final int FLAG_SECOND_FF = 129;

  static final int[] FLIP = new int[178];
  static {
    FLIP[FLAG_FIRST_FR] = FLAG_SECOND_RF;
    FLIP[FLAG_FIRST_RF] = FLAG_SECOND_FR;
    FLIP[FLAG_FIRST_RR] = FLAG_SECOND_RR;
    FLIP[FLAG_FIRST_FF] = FLAG_SECOND_FF;
  }

  static int flip(final int flag) {
    return FLIP[flag];
  }

  public static SAMRecord makeTestRecord(SAMFileHeader sfh, int flags, String sequence, String mateSequence, int pos, int matePos, String readGroup) {
    final SAMRecord sr = new SAMRecord(sfh);
    sr.setFlags(flags);
    sr.setReferenceName(sequence);
    sr.setMateReferenceName(mateSequence);
    sr.setAlignmentStart(pos);
    sr.setMateAlignmentStart(matePos);
    sr.setAttribute("RG", readGroup);
    if (pos < matePos) {
      sr.setInferredInsertSize(matePos - pos + 5);
    } else {
      sr.setInferredInsertSize(matePos - pos  - 5); //-1 * (pos - matePos + 5));
    }
    return sr;
  }

  SAMRecord orientRec(int flag) {
    final SAMFileHeader sfh = new SAMFileHeader();
    sfh.addReadGroup(new SAMReadGroupRecord("RG1"));
    return  makeTestRecord(sfh, flag, "s1", "s2", 10, 85, "RG1");
  }

  public void testOrientation() {
    assertEquals(Orientation.DU, Orientation.orientation(orientRec(FLAG_FIRST_RF), MachineOrientation.FR));
    assertEquals(Orientation.UD, Orientation.orientation(orientRec(FLAG_FIRST_FR), MachineOrientation.FR));
    assertEquals(Orientation.UU, Orientation.orientation(orientRec(FLAG_FIRST_FF), MachineOrientation.FR));
    assertEquals(Orientation.DD, Orientation.orientation(orientRec(FLAG_FIRST_RR), MachineOrientation.FR));

    assertEquals(Orientation.UD, Orientation.orientation(orientRec(FLAG_SECOND_FR), MachineOrientation.FR));
    assertEquals(Orientation.DU, Orientation.orientation(orientRec(FLAG_SECOND_RF), MachineOrientation.FR));
    assertEquals(Orientation.UU, Orientation.orientation(orientRec(FLAG_SECOND_FF), MachineOrientation.FR));
    assertEquals(Orientation.DD, Orientation.orientation(orientRec(FLAG_SECOND_RR), MachineOrientation.FR));

    assertEquals(Orientation.DU, Orientation.orientation(orientRec(FLAG_FIRST_RR), MachineOrientation.TANDEM));
    assertEquals(Orientation.UD, Orientation.orientation(orientRec(FLAG_FIRST_FF), MachineOrientation.TANDEM));
    assertEquals(Orientation.UU, Orientation.orientation(orientRec(FLAG_FIRST_FR), MachineOrientation.TANDEM));
    assertEquals(Orientation.DD, Orientation.orientation(orientRec(FLAG_FIRST_RF), MachineOrientation.TANDEM));

    assertEquals(Orientation.UD, Orientation.orientation(orientRec(FLAG_SECOND_RR), MachineOrientation.TANDEM));
    assertEquals(Orientation.DU, Orientation.orientation(orientRec(FLAG_SECOND_FF), MachineOrientation.TANDEM));
    assertEquals(Orientation.UU, Orientation.orientation(orientRec(FLAG_SECOND_RF), MachineOrientation.TANDEM));
    assertEquals(Orientation.DD, Orientation.orientation(orientRec(FLAG_SECOND_FR), MachineOrientation.TANDEM));

  }

  public void testEnum() {
    TestUtils.testEnum(Orientation.class, "[UU, UD, DU, DD]");
  }

  public void testConstructors() {
    assertEquals(+1, Orientation.UU.xDir());
    assertEquals(+1, Orientation.UU.yDir());
    assertEquals(+1, Orientation.UD.xDir());
    assertEquals(-1, Orientation.UD.yDir());
    assertEquals(-1, Orientation.DU.xDir());
    assertEquals(+1, Orientation.DU.yDir());
    assertEquals(-1, Orientation.DD.xDir());
    assertEquals(-1, Orientation.DD.yDir());
  }
  public void testXY() {
    assertEquals(+42, Orientation.UU.x(42));
    assertEquals(+42, Orientation.UU.y(42));
    assertEquals(+42, Orientation.UD.x(42));
    assertEquals(-42, Orientation.UD.y(42));
    assertEquals(-42, Orientation.DU.x(42));
    assertEquals(+42, Orientation.DU.y(42));
    assertEquals(-42, Orientation.DD.x(42));
    assertEquals(-42, Orientation.DD.y(42));

    assertEquals(+42.1, Orientation.UU.x(42.1));
    assertEquals(+42.1, Orientation.UU.y(42.1));
    assertEquals(+42.1, Orientation.UD.x(42.1));
    assertEquals(-42.1, Orientation.UD.y(42.1));
    assertEquals(-42.1, Orientation.DU.x(42.1));
    assertEquals(+42.1, Orientation.DU.y(42.1));
    assertEquals(-42.1, Orientation.DD.x(42.1));
    assertEquals(-42.1, Orientation.DD.y(42.1));
  }


  public void testFlip() {
    assertEquals(Orientation.UU, Orientation.UU.flip());
    assertEquals(Orientation.DU, Orientation.UD.flip());
    assertEquals(Orientation.UD, Orientation.DU.flip());
    assertEquals(Orientation.DD, Orientation.DD.flip());
  }

  public void testStaticMethod() {
    assertEquals(Orientation.UU, Orientation.orientation(1, 1));
    assertEquals(Orientation.UD, Orientation.orientation(1, -1));
    assertEquals(Orientation.DU, Orientation.orientation(-1, 1));
    assertEquals(Orientation.DD, Orientation.orientation(-1, -1));
    assertEquals(Orientation.UU, Orientation.orientation(0, 0));
  }

}
