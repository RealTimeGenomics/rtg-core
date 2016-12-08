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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.relation.Family;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import junit.framework.TestCase;

/**
 */
public class SegregationVcfAnnotatorTest extends TestCase {

  public void test() {
    final VcfHeader header = new VcfHeader();
    header.addSampleName("FATHER").addSampleName("MOTHER");
    final String[] childrenNames = new String[8];
    for (int i = 0; i < 4; ++i) {
      final String sonName = "SON" + (i + 1);
      final String daughterName = "DAUGHTER" + (i + 1);
      header.addSampleName(sonName).addSampleName(daughterName);
      childrenNames[i * 2] = sonName;
      childrenNames[i * 2 + 1] = daughterName;
    }
    final Family family = new Family("FATHER", "MOTHER", childrenNames);
    final SegregationVcfAnnotator annotator = new SegregationVcfAnnotator(family);
    annotator.updateHeader(header);
    annotator.updateHeader(header); // Should not be a problem doing this twice
    checkRecord("0/1\t1/1\t.\t0/1\t1/1\t1/1\t0/1\t1/1\t1/1\t1/1", annotator, "-1.808");
    checkRecord("0/1\t1/0\t0/0\t0/1\t1/1\t0/1\t0/1\t1/1\t0/0\t1/1", annotator, "-2.683");
    checkRecord("0\t0/1\t0\t0/0\t1\t0/1\t1\t0/0\t0\t0/1", annotator, "-1.297");
    checkRecord(".\t0/1\t0\t0/0\t1\t0/1\t1\t0/0\t0\t0/1", annotator, "0.000");
    checkRecord("0/1\t.\t0\t0/0\t1\t0/1\t1\t0/0\t0\t0/1", annotator, "0.000");
    checkRecord(".\t.\t0\t0/0\t1\t0/1\t1\t0/0\t0\t0/1", annotator, null);
  }

  private void checkRecord(String samples, VcfAnnotator annotator, String expected) {
    final VcfRecord rec = VcfReader.vcfLineToRecord("chr20\t1400\t.\tG\tA\t29\tPASS\t.\tGT\t" + samples);
    annotator.annotate(rec);
    annotator.annotate(rec); // Should not be a problem doing this twice
    if (expected == null) {
      assertEquals(0, rec.getInfo().size());
    } else {
      assertEquals(1, rec.getInfo().size());
      assertNotNull(rec.getInfo().get("SGP"));
      assertEquals(1, rec.getInfo().get("SGP").size());
      assertEquals(expected, rec.getInfo().get("SGP").get(0));
    }
  }

  public void testUpdateHeaderFailure() {
    final VcfHeader header = new VcfHeader();
    header.addInfoField("SGP", MetaType.INTEGER, VcfNumber.ONE, "Conflicting Header");
    header.addSampleName("FATHER").addSampleName("MOTHER").addSampleName("CHILD");
    final Family family = new Family("FATHER", "MOTHER", "CHILD");
    final SegregationVcfAnnotator annotator = new SegregationVcfAnnotator(family);
    try {
      annotator.updateHeader(header);
      fail("Allowed conflicing header line");
    } catch (NoTalkbackSlimException e) {
      //expected
      assertEquals("A VCF INFO field SGP which is incompatible is already present in the VCF header.", e.getMessage());
    }
  }

  public void testCheckHeader() {
    final VcfHeader header = new VcfHeader();
    header.addSampleName("FATHER").addSampleName("MOTHER");
    final String[] childrenNames = new String[8];
    for (int i = 0; i < 4; ++i) {
      final String sonName = "SON" + (i + 1);
      final String daughterName = "DAUGHTER" + (i + 1);
      header.addSampleName(sonName).addSampleName(daughterName);
      childrenNames[i * 2] = sonName;
      childrenNames[i * 2 + 1] = daughterName;
    }
    Family family = new Family("FATHER", "MOTHER", childrenNames);
    assertTrue(SegregationVcfAnnotator.checkHeader(header, family));
    family = new Family("FATHER", "FATHER", childrenNames);
    assertFalse(SegregationVcfAnnotator.checkHeader(header, family));
    family = new Family("FATHER", "MOTHER", "FATHER", "SON1");
    assertFalse(SegregationVcfAnnotator.checkHeader(header, family));
  }

  public void testGetSegregationScore() {
    final Code code = new CodeDiploid(4);
    ISegregationScore score = SegregationVcfAnnotator.getSegregationScore(code, -1, code.code(0, 0), false, false);
    assertTrue(score instanceof SegregationTrivial);
    score = SegregationVcfAnnotator.getSegregationScore(code, code.code(0, 0), -1, false, false);
    assertTrue(score instanceof SegregationTrivial);
    score = SegregationVcfAnnotator.getSegregationScore(code, code.code(0, 0), code.code(0, 0), true, true);
    assertTrue(score instanceof SegregationTrivial);
    score = SegregationVcfAnnotator.getSegregationScore(code, code.code(0, 0), code.code(0, 0), false, false);
    assertTrue(score instanceof SegregationScore);
    score = SegregationVcfAnnotator.getSegregationScore(code, code.code(0), code.code(0, 1), true, false);
    assertTrue(score instanceof AbstractSegregationHaploid);
    score = SegregationVcfAnnotator.getSegregationScore(code, code.code(0, 1), code.code(1), false, true);
    assertTrue(score instanceof AbstractSegregationHaploid);
  }
}
