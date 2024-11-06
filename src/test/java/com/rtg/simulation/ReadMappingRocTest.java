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
package com.rtg.simulation;

import junit.framework.TestCase;

/**
 */
public class ReadMappingRocTest extends TestCase {

  public void testROC() {
    final ReadMappingRoc roc = new ReadMappingRoc("test");
    assertEquals(0, roc.getTp(20), 0.001);
    assertEquals(0, roc.getMaxScore());

    roc.addTp(1);
    roc.addTp(1);
    roc.addTp(1);
    roc.addTp(1);
    roc.addFp(2);
    roc.addFp(2);
    roc.addFp(3);
    roc.addFp(3);

    assertEquals(3, roc.getMaxScore());

    assertEquals(0, roc.getTp(0), 0.001);
    assertEquals(0, roc.getFp(0), 0.001);

    assertEquals(4, roc.getTp(1), 0.001);
    assertEquals(0, roc.getFp(1), 0.001);

    assertEquals(0, roc.getTp(2), 0.001);
    assertEquals(2, roc.getFp(2), 0.001);

    assertEquals(0, roc.getTp(3), 0.001);
    assertEquals(2, roc.getFp(3), 0.001);

    roc.addFp(3, 0.25);
    assertEquals(2.25, roc.getFp(3), 0.001);

    //System.err.println(roc.getDistribution());
    //System.err.println(roc.getRoc());
  }

}
