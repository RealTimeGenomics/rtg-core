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
import com.rtg.util.TestUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;

/**
 */
public class KappaImplementationTest extends AbstractKappaTest {

  public void testToString() throws InvalidParamsException, IOException {
    final AbstractKappa kappa;
    kappa = getKappa();
    //System.out.println(kappa.toString());
    final String exp = FileHelper.resourceToString("com/rtg/variant/bayes/complex/resources/kappatest0.txt");
    assertTrue(TestUtils.sameLines(exp, kappa.toString(), false));
  }

  public void testBounds() {
    final double iota = 0.02;
    final double delta = 0.03;
    final AbstractMachineErrorParams pp = MachineErrorParams.builder()
        .errorInsEventRate(iota)
        .errorInsDistribution(new double[] {1.0, 0.0})
        .errorDelEventRate(delta)
        .errorDelDistribution(new double[] {1.0, 0.0})
        .create();
    final KappaImplementation ip = new KappaImplementation(pp, 0.0);
    try {
      ip.kappa(0, -1);
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("-1", e.getMessage());
    }
    try {
      ip.kappa(-1, 0);
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("-1", e.getMessage());
    }
  }

  public void testCorrectionKappa() {
    final double[] corr = KappaImplementation.createCumuPi(new double[] {0.1, 0.7, 0.2});
    assertEquals(3, corr.length);
    assertEquals(1.0, corr[0], 0.000001);
    assertEquals(0.9, corr[1], 0.000001);
    assertEquals(0.2, corr[2], 0.000001);
  }
}
