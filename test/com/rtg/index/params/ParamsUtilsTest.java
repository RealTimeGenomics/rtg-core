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
package com.rtg.index.params;

import junit.framework.TestCase;

/**
 */
public class ParamsUtilsTest extends TestCase {

  public final void testMemToString() {
    assertEquals("\tMemory\tfoo\t42" + com.rtg.util.StringUtils.LS, ParamsUtils.memToString("foo", 42));
    try {
      ParamsUtils.memToString("foo bar", 42);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
  }

  public final void testMemToString2() {
    assertEquals("\tMemory\tfoo\t42\t21\t-1" + com.rtg.util.StringUtils.LS, ParamsUtils.memToString("foo", 42, 21, -1));
  }

}

