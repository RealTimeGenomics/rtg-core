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

package com.rtg.segregation;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class SearchTypeTest extends TestCase {

  public void test() {
    TestUtils.testEnum(SearchType.class, "[XO, New, Error, OK]");
    assertEquals("BX", SearchType.XO.code());
    assertEquals("BN", SearchType.New.code());
    assertEquals("BL", SearchType.OK.code());
    assertEquals("BE", SearchType.Error.code());

    assertEquals("X", SearchType.XO.bedCode());
    assertEquals("N", SearchType.New.bedCode());
    assertEquals("N", SearchType.OK.bedCode());
    assertEquals(null, SearchType.Error.bedCode());
  }
}
