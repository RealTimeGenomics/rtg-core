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
package com.rtg.index.hash.ngs;

import com.rtg.index.hash.ngs.instances.CGMaska1b1;
import com.rtg.index.hash.ngs.instances.MaskL36w18s3e1;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;

import junit.framework.TestCase;

/**
 *
 *
 */
public class FactoryUtilTest extends TestCase {

  /**
   */
  public FactoryUtilTest(final String name) {
    super(name);
  }

  public void testFactory() {
    final HashFunctionFactory factory = FactoryUtil.hashFunction("MaskL36w18s3e1");

    assertEquals(MaskL36w18s3e1.FACTORY.hashBits(), factory.hashBits());
    assertEquals(MaskL36w18s3e1.FACTORY.numberWindows(), factory.numberWindows());
    assertEquals(MaskL36w18s3e1.FACTORY.hashBits(), factory.hashBits());
  }

  public void testErrorMessage() {
    Diagnostic.setLogStream();
    try {
      FactoryUtil.hashFunction("foo");
      fail();
    } catch (final SlimException e) {
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testCGMask() {
    final HashFunctionFactory factory = FactoryUtil.hashFunction("CGMaska1b1");

    assertEquals(CGMaska1b1.FACTORY.hashBits(), factory.hashBits());
    assertEquals(CGMaska1b1.FACTORY.numberWindows(), factory.numberWindows());
    assertEquals(CGMaska1b1.FACTORY.hashBits(), factory.hashBits());
  }

  public void testFactories() {
    final String[] masks = {"MaskL30w12s3e1", "MaskL35w15s2e1", "MaskL36w18s3e1",
        "MaskL51w12s3e3", "MaskL51w15s3e2", "MaskL51w17s2e2", "MaskL51w18s2e1", "SplitL36w12s2e2",
        "SplitL36w12s3e2", "SplitL36w12s3e3", "SplitL36w18s1e1", "SplitL36w18s2e1", "SplitL4w2s1e1",
        "SplitL4w2s1e1b", "SplitL4w4s0e0", "CGMaska0", "CGMaska1b1", "CGMaska15b1"};
    for (String mask : masks) {
      assertTrue(FactoryUtil.checkMaskExists(mask));
    }
  }

}
