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

package com.rtg.variant.sv;

import junit.framework.TestCase;

/**
 */
public class DistributionUtilsTest extends TestCase {

  public void testDistribution1() {
    final Distribution di = DistributionUtils.distribution1(-5, 5, 1.0, -2.0, 2.0);
    di.globalIntegrity();
    //System.err.println(di.dump());
    assertEquals(0.9332, di.get(-5), 0.0001);
    assertEquals(0.5, di.get(-2), 0.0001);
    assertEquals(0.0013, di.get(4), 0.0001);
  }

  public void testDistribution12() {
    final Distribution di = DistributionUtils.distribution1(-5, 5, 2.0, -2.0, 2.0);
    di.globalIntegrity();
    //System.err.println(di.dump());
    assertEquals(0.9332 * 2, di.get(-5), 0.0001);
    assertEquals(0.5 * 2.0, di.get(-2), 0.0001);
    assertEquals(0.0013 * 2, di.get(4), 0.0001);
  }

  public void testDistribution2() {
    final Distribution di = DistributionUtils.distribution2(-5, 5, 1.0, -2.0, 2.0);
    di.globalIntegrity();
    //System.err.println(di.dump());
    assertEquals(0.0668, di.get(-5), 0.0001);
    assertEquals(0.5, di.get(-2), 0.0001);
    assertEquals(0.9987, di.get(4), 0.0001);
  }

  public void testDistribution22() {
    final Distribution di = DistributionUtils.distribution2(-5, 5, 2.0, -2.0, 2.0);
    di.globalIntegrity();
    //System.err.println(di.dump());
    assertEquals(0.0668 * 2.0, di.get(-5), 0.0001);
    assertEquals(0.5 * 2.0, di.get(-2), 0.0001);
    assertEquals(0.9987 * 2.0, di.get(4), 0.0001);
  }

  public void testAdd1() {
    final Distribution su = DistributionUtils.add(new DistributionArray(-2, new double[] {0.9, 0.8, 0.0, 0.0}), 1.0);
    su.globalIntegrity();
    assertEquals(1.9, su.get(-2), 0.0000001);
    assertEquals(1.8, su.get(-1), 0.0000001);
    assertEquals(1.0, su.get(0), 0.0000001);
    assertEquals(1.0, su.get(1), 0.0000001);
  }

  public void testAdd2() {
    final Distribution su = DistributionUtils.add(new DistributionArray(-2, new double[] {0.9, 0.8, 0.0, 0.0}), new DistributionArray(-2, new double[] {0.9, 0.8, 0.0, 0.0}));
    su.globalIntegrity();
    assertEquals(1.8, su.get(-2), 0.0000001);
    assertEquals(1.6, su.get(-1), 0.0000001);
    assertEquals(0.0, su.get(0), 0.0000001);
    assertEquals(0.0, su.get(1), 0.0000001);
  }

  public void testMultiply() {
    final Distribution su = DistributionUtils.multiply(new DistributionArray(-2, new double[] {0.9, 0.8, 1.0, 1.0}), new DistributionArray(-2, new double[] {0.9, 0.8, 0.0, 0.0}));
    su.globalIntegrity();
    assertEquals(0.81, su.get(-2), 0.0000001);
    assertEquals(0.64, su.get(-1), 0.0000001);
    assertEquals(0.0, su.get(0), 0.0000001);
    assertEquals(0.0, su.get(1), 0.0000001);
  }

  public void testSubtract() {
    final Distribution su = DistributionUtils.subtract(1.0, new DistributionArray(-2, new double[] {0.9, 0.8, 0.0, 0.0}));
    su.globalIntegrity();
    assertEquals(0.1, su.get(-2), 0.0000001);
    assertEquals(0.2, su.get(-1), 0.0000001);
    assertEquals(1.0, su.get(0), 0.0000001);
    assertEquals(1.0, su.get(1), 0.0000001);
  }

  public void testSubtractEdge() {
    final Distribution su = DistributionUtils.subtract(1.0, new DistributionArray(-2, new double[] {1.0000001, 0.8, 0.0, 0.0}));
    su.globalIntegrity();
    assertEquals(0.0, su.get(-2), 0.0000001);
    assertEquals(0.2, su.get(-1), 0.0000001);
    assertEquals(1.0, su.get(0), 0.0000001);
    assertEquals(1.0, su.get(1), 0.0000001);
  }

  public void testSubtractBad() {
    try {
      DistributionUtils.subtract(1.0, new DistributionArray(-2, new double[] {1.01, 0.8, 0.0, 0.0}));
      fail();
    } catch (final RuntimeException e) {
    }
  }
}
