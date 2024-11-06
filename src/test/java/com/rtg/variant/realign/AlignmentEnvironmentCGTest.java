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

package com.rtg.variant.realign;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class AlignmentEnvironmentCGTest extends TestCase {

  public void test() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/2");
    rec.setReadPairedFlag(true);
    final VariantParams params = VariantParams.builder().create();
    final AbstractMachineErrorParams me = MachineErrorParams.builder().machine(MachineType.COMPLETE_GENOMICS).create();
    final AlignmentEnvironment se = new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, null, me.machineType());

    assertFalse(!se.isInverted());
    final String exp = ""
      + "AlignmentEnvironment read=GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT "
      + "quality=["
      + "0.0200, 0.0398, 0.0126, 0.0100, 0.0200, 0.0200, 0.0200, 0.0200, 1.0000, 0.0200, "
      + "0.0100, 0.0501, 0.0794, 0.0316, 0.0794, 0.0050, 0.0079, 0.0501, 0.0251, 0.0100, "
      + "0.0200, 0.0100, 0.0200, 0.0100, 0.0631, 0.0316, 0.0126, 0.0050, 0.0251, 0.0079, "
      + "0.0158, 0.0100, 0.0200, 0.0063, 0.0398"
      + "] start=-1";
    assertEquals(exp, se.toString());
  }

  public void testDefaultQuality() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    rec.setReadPairedFlag(true);
    final VariantParams params = VariantParams.builder().defaultQuality(10).create();
    final AlignmentEnvironment se = new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, null, MachineErrorParams.builder().machine(MachineType.COMPLETE_GENOMICS).create().machineType());

    assertFalse(!se.isInverted());
    final String exp = ""
      + "AlignmentEnvironment read=GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT "
      + "quality=["
      + "0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, "
      + "0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, "
      + "0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, "
      + "0.1000, 0.1000, 0.1000, 0.1000, 0.1000"
      + "] start=-1";
    assertEquals(exp, se.toString());
  }

  public void testBad() {
    Diagnostic.setLogStream();
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/2");
    final VariantParams params = VariantParams.builder().create();
    try {
      new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, null, MachineErrorParams.builder().machine(MachineType.COMPLETE_GENOMICS).create().machineType());
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertTrue(e.getMessage(), e.getMessage().startsWith("Invalid CG alignment."));
    }
  }
}
