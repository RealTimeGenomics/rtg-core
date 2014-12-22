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
package com.rtg.reader;

import junit.framework.TestCase;

/**
 * Test class for ArrayPrereadNames
 */
public class ArrayPrereadNamesTest extends TestCase {

  public void testPrereadNames() {
    final ArrayPrereadNames pn = new ArrayPrereadNames(new String[] {"blah", "fah", "rah"});
    assertEquals("fah", pn.name(1));
    assertEquals(3, pn.length());
  }

}
