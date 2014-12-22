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

import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;

/**
 */
public class KappaMakeTester extends AbstractKappaTest {

  @Override
  protected AbstractKappa getKappa() {
    final AbstractMachineErrorParams pp = MachineErrorParams.builder()
    .errorDelEventRate(0.03)
    .errorDelDistribution(new double[] {0.8, 0.2})
    .errorInsEventRate(0.02)
    .errorInsDistribution(new double[] {0.75, 0.25})
    .create();
    final AbstractKappa kappa = AbstractKappa.makeKappa(pp.errorInsEventRate(), pp.errorInsDistribution(), pp.errorDelEventRate(), pp.errorDelDistribution(), 0.2);
    kappa.globalIntegrity();
    return kappa;
  }

}
