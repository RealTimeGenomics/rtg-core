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
package com.rtg.variant.sv;

import com.rtg.util.machine.MachineType;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class ReadGroupStateTest extends TestCase {


  public void testState() {
    final ReadGroupStats stats = new ReadGroupStats("1", 3);

    final ReadGroupState rgs = new ReadGroupState(stats, MachineType.COMPLETE_GENOMICS, null);

    assertEquals(0, rgs.notPaired().length());
    assertEquals(0, rgs.unmatedLeftArm(false).length());
    assertEquals(0, rgs.unmatedRightArm(false).length());
    assertEquals(0, rgs.properLeftArm(false).length());
    assertEquals(0, rgs.properRightArm(false).length());
    assertEquals(0, rgs.discordantLeftArm(false).length());
    assertEquals(0, rgs.discordantRightArm(false).length());
    assertEquals(0, rgs.unique().length());
    assertEquals(0, rgs.ambiguous().length());
    assertEquals(rgs.unmatedLeftArm(false), rgs.unmatedRightArm(true));
    assertEquals(rgs.unmatedLeftArm(true), rgs.unmatedRightArm(false));
    assertEquals(rgs.properLeftArm(false), rgs.properRightArm(true));
    assertEquals(rgs.properLeftArm(true), rgs.properRightArm(false));
    assertEquals(rgs.discordantLeftArm(false), rgs.discordantRightArm(true));
    assertEquals(rgs.discordantLeftArm(true), rgs.discordantRightArm(false));

    assertEquals(0, rgs.lo());
    assertEquals(1, rgs.hi());

    assertEquals(stats, rgs.stats());
    assertEquals(MachineType.COMPLETE_GENOMICS, rgs.machine());
    assertEquals(MachineType.COMPLETE_GENOMICS.orientation(), rgs.orientation());

    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    r.setAlignmentStart(5);
    r.setCigarString("2=");

    rgs.reset(5, 5);
    assertTrue(rgs.update(r));

    r.setFlags(153);
    assertTrue(rgs.update(r));

    r.setFlags(131);
    assertTrue(rgs.update(r));

    assertEquals(5, rgs.notPaired().length());
    assertEquals(5, rgs.unmatedLeftArm(false).length());
    assertEquals(5, rgs.unmatedRightArm(false).length());
    assertEquals(5, rgs.properLeftArm(false).length());
    assertEquals(5, rgs.properRightArm(false).length());
    assertEquals(5, rgs.discordantLeftArm(false).length());
    assertEquals(5, rgs.discordantRightArm(false).length());
    assertEquals(5, rgs.unique().length());
    assertEquals(5, rgs.ambiguous().length());

    assertEquals(0, rgs.lo());
    assertEquals(1, rgs.hi());

  }

}
