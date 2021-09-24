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
package com.rtg.variant;

import junit.framework.TestCase;

/**
 * Test class for {@link NoQualityException}
 */
public class NoQualityExceptionTest extends TestCase {

  /**
   */
  public NoQualityExceptionTest(final String name) {
    super(name);
  }

  public void test() {
    final NoQualityException e = new NoQualityException();
    assertNotNull(e);
    //Don't check message, it changes on platforms and CLRs.
    //String expected = null;
    //assertEquals(expected, e.getMessage());
  }
}
