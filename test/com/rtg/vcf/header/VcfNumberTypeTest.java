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
package com.rtg.vcf.header;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class VcfNumberTypeTest extends TestCase {

  public void testSomeMethod() {
    TestUtils.testEnum(VcfNumberType.class, "[INTEGER, A, R, G, .]");
    assertTrue(VcfNumberType.ALTS.isFixedString());
    assertTrue(VcfNumberType.GENOTYPES.isFixedString());
    assertTrue(VcfNumberType.UNKNOWN.isFixedString());
    assertFalse(VcfNumberType.INTEGER.isFixedString());
  }
}
