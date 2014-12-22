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

/**
 */
public class KappaMemoTest extends AbstractKappaTest {

  @Override
  protected AbstractKappa getKappa() {
    return new KappaMemo(super.getKappa());
  }

  public void testToString() throws InvalidParamsException, IOException {
    final AbstractKappa kappa;
    kappa = getKappa();
    //System.out.println(kappa.toString());
    final String exp = FileHelper.resourceToString("com/rtg/variant/bayes/complex/resources/kappamemotest0.txt");
    assertTrue(TestUtils.sameLines(exp, kappa.toString(), false));
  }

}
