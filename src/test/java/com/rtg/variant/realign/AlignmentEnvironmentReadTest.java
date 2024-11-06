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

import com.rtg.util.machine.MachineType;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class AlignmentEnvironmentReadTest extends TestCase {

  public void test() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString("NACGT");
    sam.setCigarString("5=");
    sam.setBaseQualityString("!!0AB");
    sam.setAlignmentStart(42);
    final VariantParams params = VariantParams.builder().create();
    final AbstractAlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineType.ILLUMINA_PE);
    se.integrity();
    assertEquals("AlignmentEnvironment read=NACGT quality=[1.0000, 1.0000, 0.0316, 0.0006, 0.0005] start=41", se.toString());
    assertEquals(0, se.base(0));
    assertEquals(1, se.base(1));
    assertEquals(2, se.base(2));
    assertEquals(3, se.base(3));
    assertEquals(4, se.base(4));

    assertEquals(1.0, se.quality(0));
    assertEquals(1.0, se.quality(1));
    assertEquals(0.0316, se.quality(2), 0.0001);
    assertEquals(0.0006, se.quality(3), 0.0001);
    assertEquals(0.0005, se.quality(4), 0.0001);

    assertEquals(41, se.start());

    assertEquals(5, se.subsequenceLength());

    assertTrue(!se.isInverted());
  }

  public void testBug() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString("CCCCCTAGGGGG");
    sam.setCigarString("5=2I5=");
    sam.setAlignmentStart(1);
    final VariantParams params = VariantParams.builder().create();
    final AbstractAlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineType.ILLUMINA_PE);
    se.integrity();
    assertEquals("AlignmentEnvironment read=CCCCCTAGGGGG quality=[0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100] start=0", se.toString());
    assertEquals(2, se.base(0));
    assertEquals(2, se.base(4));
    assertEquals(4, se.base(5));
    assertEquals(1, se.base(6));
    assertEquals(3, se.base(7));
    assertEquals(3, se.base(11));

    assertEquals(0.0100, se.quality(0), 0.0000001);
    assertEquals(0.0100, se.quality(11), 0.0000001);

    assertEquals(0, se.start());

    assertEquals(12, se.subsequenceLength());
  }

}
