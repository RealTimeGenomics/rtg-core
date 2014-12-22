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
package com.rtg.util.array.zeroindex;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;


/**
 */
public class ZeroCreateTest extends TestCase {

  public void testOK() {
    final ZeroIndex index = ZeroCreate.createIndex(10L);
    assertEquals("Index [10]" + LS, index.toString());
  }

}
