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
    final String[] masks = {"CGMaska0", "CGMaska1b1", "CGMaska15b1"};
    for (String mask : masks) {
      assertTrue(FactoryUtil.checkMaskExists(mask));
    }
  }

}
