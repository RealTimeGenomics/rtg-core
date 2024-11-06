/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
