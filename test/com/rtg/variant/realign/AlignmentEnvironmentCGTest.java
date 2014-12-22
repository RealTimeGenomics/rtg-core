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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

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
    final AlignmentEnvironment se = new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, me, null);

    assertFalse(se.cgOverlapOnLeft());
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
    final AlignmentEnvironment se = new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, MachineErrorParams.builder().machine(MachineType.COMPLETE_GENOMICS).create(), null);

    assertFalse(se.cgOverlapOnLeft());
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
      new AlignmentEnvironmentCG(new VariantAlignmentRecord(rec), params, MachineErrorParams.builder().machine(MachineType.COMPLETE_GENOMICS).create(), null);
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertTrue(e.getMessage(), e.getMessage().startsWith("Invalid CG alignment record="));
    }
  }
}
