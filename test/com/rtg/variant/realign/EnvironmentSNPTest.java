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

package com.rtg.variant.realign;

import junit.framework.TestCase;

/**
 */
public class EnvironmentSNPTest extends TestCase {

  public void test1() {
    final Environment env0 = new EnvironmentImplementation(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        3, //start
        new byte[] {1, 2, 3, 4, 0}, //read
        new double[]{0.1, 0.001, 0.002, 0.003, 0.25} //quality
    );
    final Environment env = new EnvironmentSNP(env0, 2, (byte) 4);
    assertEquals(7, env.templateLength());
    assertEquals(0, env.template(-4));
    assertEquals(1, env.template(-3));
    assertEquals(4, env.template(0));
    assertEquals(0, env.template(1));
    assertEquals(4, env.template(2));  // this one has changed!
    assertEquals(2, env.template(3));
    assertEquals(0, env.template(4));
  }
}
