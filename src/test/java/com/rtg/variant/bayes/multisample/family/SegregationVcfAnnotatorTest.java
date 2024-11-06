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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.relation.Family;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfReaderTest;
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
    final VcfRecord rec = VcfReaderTest.vcfLineToRecord("chr20\t1400\t.\tG\tA\t29\tPASS\t.\tGT\t" + samples);
    annotator.annotate(rec);
    annotator.annotate(rec); // Should not be a problem doing this twice
    if (expected == null) {
      assertEquals(0, rec.getInfo().size());
    } else {
      assertEquals(1, rec.getInfo().size());
      assertEquals(expected, rec.getInfo("SGP"));
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
