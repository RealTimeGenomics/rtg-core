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
package com.rtg.variant.bayes.complex;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class IndelUtilsTest extends TestCase {

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#flatten(double[], double)}.
   */
  public final void testFlatten1() {
    final double[] f = IndelUtils.flatten(new double[] {0.0, 0.5, 0.25, 0.25}, 0.3);
    assertEquals("[0.700, 0.150, 0.075, 0.075]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#flatten(double[], double)}.
   */
  public final void testFlatten2() {
    final double[] f = IndelUtils.flatten(new double[] {0.0, 1.0}, 0.3);
    assertEquals("[0.700, 0.300]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#unFlatten(double[])}.
   */
  public final void testUnFlatten1() {
    final double[] f = IndelUtils.unFlatten(new double[] {0.2, 0.2, 0.2, 0.2, 0.2});
    assertEquals("[0.000, 0.250, 0.250, 0.250, 0.250]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#unFlatten(double[])}.
   */
  public final void testUnFlatten2() {
    final double[] f = IndelUtils.unFlatten(new double[] {0.5, 0.5});
    assertEquals("[0.000, 1.000]", Utils.realFormat(f, 3));
  }
}
