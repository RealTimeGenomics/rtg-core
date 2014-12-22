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
package com.rtg.variant.bayes.complex;

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractKappaTest extends TestCase {

  public void testKappaPriors() throws InvalidParamsException, IOException {
    final AbstractKappa kappa;
    kappa = getKappa();

    assertEquals(0.0012, kappa.pi(-3), 0.000001);
    assertEquals(0.0060, kappa.pi(-2), 0.000001);
    assertEquals(0.0240, kappa.pi(-1), 0.000001);
    assertEquals(0.9500, kappa.pi(0), 0.000001);
    assertEquals(0.0150, kappa.pi(1), 0.000001);
    assertEquals(0.0050, kappa.pi(2), 0.000001);
    assertEquals(0.0010, kappa.pi(3), 0.000001);
    assertEquals(0.0002, kappa.pi(4), 0.000001);

    assertEquals(1.00245, kappa.piSum(-3), 0.000001);
    assertEquals(1.00125, kappa.piSum(-2), 0.000001);
    assertEquals(0.99525, kappa.piSum(-1), 0.000001);
    assertEquals(0.97125, kappa.piSum(0), 0.000001);
    assertEquals(0.02125, kappa.piSum(1), 0.000001);
    assertEquals(0.00625, kappa.piSum(2), 0.000001);
    assertEquals(0.00125, kappa.piSum(3), 0.000001);
    assertEquals(0.00025, kappa.piSum(4), 0.000001);

    check(kappa, 0, 0);
    assertEquals(0.978120978, kappa.kappa(0, 0), 0.000001);
    assertEquals(0.001197067, kappa.kappa(0, 3), 0.000001);
    assertEquals(0.001029601, kappa.kappa(3, 0), 0.000001);
    assertEquals(0.000205878, kappa.kappa(4, 0), 0.000001);
    assertEquals(0.947678188, kappa.kappa(3, 3), 0.000001);
  }

  protected AbstractKappa getKappa() throws InvalidParamsException {
    final AbstractKappa kappa;
    final AbstractMachineErrorParams pp = MachineErrorParams.builder()
    .errorDelEventRate(0.03)
    .errorDelDistribution(new double[] {0.8, 0.2})
    .errorInsEventRate(0.02)
    .errorInsDistribution(new double[] {0.75, 0.25})
    .create();
    kappa = new KappaImplementation(pp, 0.2);
    kappa.globalIntegrity();
    return kappa;
  }

  private void check(final AbstractKappa kappa, final int m, final int l) {
    assertEquals(kappa.kappa(m, l), kappa.pi(m - l), kappa.piSum(-l));
  }
}
