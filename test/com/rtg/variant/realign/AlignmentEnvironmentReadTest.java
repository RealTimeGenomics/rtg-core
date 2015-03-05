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

package com.rtg.variant.realign;

import com.rtg.variant.MachineErrorParams;
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
    final AbstractAlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineErrorParams.builder().create());
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

    assertTrue(se.cgOverlapOnLeft());
  }

  public void testBug() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString("CCCCCTAGGGGG");
    sam.setCigarString("5=2I5=");
    sam.setAlignmentStart(1);
    final VariantParams params = VariantParams.builder().create();
    final AbstractAlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineErrorParams.builder().create());
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
