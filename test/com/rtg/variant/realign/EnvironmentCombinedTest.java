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

import com.rtg.mode.DNA;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.complex.ComplexTemplate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class EnvironmentCombinedTest extends TestCase {

  public void test0() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString("NACGT");
    sam.setCigarString("5=");
    sam.setBaseQualityString("!!0AB");
    sam.setAlignmentStart(42);
    //System.err.println("start=" + sam.getAlignmentStart() + " end=" + sam.getAlignmentEnd());
    final VariantParams params = VariantParams.builder().create();
    final AlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineType.ILLUMINA_PE);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    //final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), 0, cot, new byte[] {});
    final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(0, 0, cot, new byte[] {});
    final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), 0, temEnv);
    env.integrity();
    assertEquals(0, env.maxShift());
    assertEquals(5, env.readLength());

    assertEquals(0, env.read(0));
    assertEquals(4, env.read(4));

    assertEquals(1.0, env.quality(0));
    assertEquals(0.0005, env.quality(4), 0.0001);

    assertEquals(41, env.absoluteTemplatePosition(0));
    assertEquals(45, env.absoluteTemplatePosition(4));
    assertEquals(0, env.templateLength());

    assertEquals(0, env.template(-1));
    assertEquals(0, env.template(0));
    assertEquals(0, env.template(1));
  }

  public void test1() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString("NACGT");
    sam.setCigarString("5=");
    sam.setBaseQualityString("!!0AB");
    sam.setAlignmentStart(1);
    final VariantParams params = VariantParams.builder().create();
    final AlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(sam), params, MachineType.ILLUMINA_PE);
    final byte[] template = {1, 1, 2, 2, 3, 3, 1, 1};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 2, 4);
    //final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), 3, cot, new byte[] {4});
    final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(0, template.length, cot, new byte[] {4});
    final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), 3, temEnv);
    //System.err.println(env);
    env.integrity();
    assertEquals(3, env.maxShift());
    assertEquals(5, env.readLength());

    assertEquals(0, env.read(0));
    assertEquals(4, env.read(4));

    assertEquals(1.0, env.quality(0));
    assertEquals(0.0005, env.quality(4), 0.0001);

    assertEquals(0, env.absoluteTemplatePosition(0));
    assertEquals(8, env.templateLength());

    assertEquals(0, env.template(-1));
    assertEquals(1, env.template(0));
    assertEquals(1, env.template(1));
    assertEquals(4, env.template(2));
    assertEquals(3, env.template(3));
    assertEquals(3, env.template(4));
    assertEquals(1, env.template(5));
    assertEquals(1, env.template(6));
    assertEquals(0, env.template(7));
  }

  private VariantParams getSpecifiedParams() {
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeIndelEventFraction(0.7)
    .genomeIndelEventRate(0.025)
    .genomeIndelDistribution(new double[] {0.5, 0.25, 0.25})
    .create();
    return new VariantParamsBuilder().genomePriors(vpb).defaultQuality(10).create();
  }

  public void testBug1() {
    final byte[] template = DNA.stringDNAtoByte("CCCCCGGGGG");
    final ComplexTemplate cot = new ComplexTemplate(template, "", 5, 5);
    final VariantParams vp = getSpecifiedParams();
    final SAMRecord samTA = new SAMRecord(new SAMFileHeader());
    samTA.setAlignmentStart(1);
    samTA.setReadString("CCCCCTAGGGGG");
    samTA.setCigarString("5=2I5=");
    final AlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(samTA), vp, MachineType.ILLUMINA_PE);

    //final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), Utils.calculateDefaultMaxShift(samTA.getReadLength()), cot, DNA.stringDNAtoByte("TA"));
    final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(0, template.length, cot, DNA.stringDNAtoByte("TA"));
    final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), MaxShiftUtils.calculateDefaultMaxShift(samTA.getReadLength()), temEnv);

    //System.err.println(env);
    final byte[] exp = {2, 2, 2, 2, 2, 4, 1, 3, 3, 3, 3, 3};
    for (int i = 0; i < 12; ++i) {
      assertEquals(exp[i], env.template(i));
      //System.err.println("[" + i + "]" + env.template(i));
    }
  }

  public void testBug2() {
    final byte[] template = DNA.stringDNAtoByte("CCCCCGGGGG");
    final ComplexTemplate cot = new ComplexTemplate(template, "", 5, 5);
    final VariantParams vp = getSpecifiedParams();
    final SAMRecord samTA = new SAMRecord(new SAMFileHeader());
    samTA.setAlignmentStart(3);
    samTA.setReadString("CCCAGGGGG");
    samTA.setCigarString("3=1I5=");
    final AlignmentEnvironment se = new AlignmentEnvironmentRead(new VariantAlignmentRecord(samTA), vp, MachineType.ILLUMINA_PE);

    //final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), Utils.calculateDefaultMaxShift(samTA.getReadLength()), cot, DNA.stringDNAtoByte("A"));
    final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(2, template.length, cot, DNA.stringDNAtoByte("A"));
    final EnvironmentCombined env = new EnvironmentCombined(se, se.start(), MaxShiftUtils.calculateDefaultMaxShift(samTA.getReadLength()), temEnv);

    //System.err.println(env);
    final byte[] exp = {2, 2, 2, 2, 2, 1, 3, 3, 3, 3, 3};
    for (int i = -2; i < 9; ++i) {
      assertEquals("" + i, exp[i + 2], env.template(i));
      //System.err.println("[" + i + "]" + env.template(i));
    }
  }
}
