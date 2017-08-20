/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.segment;

import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class CnaTypeTest extends TestCase {

  public void test() {
    TestUtils.testEnum(CnaType.class, "[DEL, DUP, NONE]");
    assertEquals("SVTYPE", VcfUtils.INFO_SVTYPE);
    assertEquals("END", VcfUtils.INFO_END);
  }

  public void testFromVcf() {
    final VcfRecord rec = new VcfRecord("seq1", 42, "A");
    assertEquals(CnaType.NONE, CnaType.valueOf(rec));
    rec.addInfo("SVTYPE", "DUP");
    assertEquals(CnaType.DUP, CnaType.valueOf(rec));
    final VcfRecord rec2 = new VcfRecord("seq1", 42, "A");
    rec2.addInfo("SVTYPE", "DEL");
    assertEquals(CnaType.DEL, CnaType.valueOf(rec2));
    final VcfRecord rec3 = new VcfRecord("seq1", 42, "A");
    rec3.addInfo("SVTYPE", "FOO");
    assertEquals(CnaType.NONE, CnaType.valueOf(rec3));
  }
}
