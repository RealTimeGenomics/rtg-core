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
