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
package com.rtg.position.output;

import junit.framework.TestCase;

/**
 */
public class UtilityTest extends TestCase {

  public void testEFormat() {
    assertEquals("0.0", Utility.eFormat(0.0));
    assertEquals("1.0", Utility.eFormat(1.0));
    assertEquals("0.12", Utility.eFormat(0.123));
    assertEquals("0.001", Utility.eFormat(0.001));
    assertEquals("1e-07", Utility.eFormat(0.000000123));
    assertEquals("20", Utility.eFormat(20.123456));
    assertEquals("9.0", Utility.eFormat(8.9567));
    assertEquals("0.90", Utility.eFormat(0.89567));
    assertEquals("0.100", Utility.eFormat(0.099678));
    assertEquals("e-110", Utility.eFormat(1.23e-110));
    assertEquals("0.0", Utility.eFormat(1.23e-200));
  }

}
