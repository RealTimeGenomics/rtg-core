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
