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
package com.rtg.variant.sv;

import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

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
