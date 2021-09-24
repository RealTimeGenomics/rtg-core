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
package com.rtg.ml;

import junit.framework.TestCase;

/**
 */
public class DatasetTest extends TestCase {

  public DatasetTest(final String name) {
    super(name);
  }

  public void test() {
    final Attribute a = new Attribute("higgs", MlDataType.DOUBLE);
    final Dataset d = new Dataset(a);
    d.addInstance(new Instance(new double[] {Math.PI}, true));
    d.addInstance(new Instance(new double[] {0.0}, true));
    d.addInstance(new Instance(new double[]{4.2}, false));
    assertEquals(a, d.getAttributes()[0]);
    assertEquals(2.0, d.totalPositiveWeight());
    assertEquals(1.0, d.totalNegativeWeight());
    assertEquals(2, d.totalPositives());
    assertEquals(1, d.totalNegatives());
    d.addInstance(new Instance(new double[]{4.2}, false, 5.0));
    assertEquals(6.0, d.totalNegativeWeight());
    assertEquals(2, d.totalPositives());
    assertEquals(2, d.totalNegatives());
    assertEquals(Math.PI, d.getInstances().get(0).instance()[0]);
    assertEquals(0.0, d.getInstances().get(1).instance()[0]);
    assertEquals(4.2, d.getInstances().get(2).instance()[0]);
    assertEquals(5.0, d.getInstances().get(3).weight());
    d.reweight();
    assertEquals(2, d.totalPositives());
    assertEquals(2, d.totalNegatives());
    assertEquals(2.0, d.totalPositiveWeight());
    assertEquals(2.0, d.totalNegativeWeight());
  }

}
