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
