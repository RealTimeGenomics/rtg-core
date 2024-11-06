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

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class EnvironmentHomopolymerTest extends TestCase {


  protected EnvironmentHomopolymer createEnv(int maxShift, byte[] template, int start, byte[] read, double[] quality) {
    return new EnvironmentHomopolymer(new EnvironmentImplementation(maxShift, template, start, read, quality));
  }

  public void test1() {
    final EnvironmentHomopolymer env = createEnv(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 1, 1, 2}, //template
        3, //start
        new byte[] {1, 2, 2, 4, 0}, //read
        new double[]{0.1, 0.001, 0.002, 0.003, 0.25} //quality
        );
    Exam.integrity(env);
    Exam.globalIntegrity(env);
    assertEquals(2, env.maxShift());
    final String exp = ""
        + "Environment" + LS
        + "Template 3..7" + LS
        + "[-2]C 1 1" + LS
        + "[-1]G 1 1" + LS
        + "[0]T 1 1" + LS
        + "[1]A 2 0" + LS
        + "[2]A 0 2" + LS
        + "[3]C 1 1" + LS
        + "Read 0..5" + LS
        + "[0]A 0.100 1 1" + LS
        + "[1]C 0.001 2 0" + LS
        + "[2]C 0.002 0 2" + LS
        + "[3]T 0.003 1 1" + LS
        + "[4]N 0.250 1 1" + LS
        ;
    assertEquals(exp, env.toString());
    assertEquals(5, env.readLength());
    assertEquals(3, env.absoluteTemplatePosition(0));
    assertEquals(7, env.absoluteTemplatePosition(4));

    assertEquals(0.1, env.quality(0));
    assertEquals(0.25, env.quality(4));

    assertEquals(1, env.read(0));
    assertEquals(0, env.read(4));

    assertEquals(7, env.templateLength());
    assertEquals(0, env.template(-4));
    assertEquals(1, env.template(-3));
    assertEquals(4, env.template(0));
    assertEquals(2, env.template(3));
    assertEquals(0, env.template(4));

    //Homopolymer stuff
    assertEquals(2, env.templateStart(1));
    assertEquals(0, env.templateEnd(1));
    assertEquals(0, env.templateStart(2));
    assertEquals(2, env.templateEnd(2));

    assertEquals(2, env.readStart(1));
    assertEquals(0, env.readEnd(1));
    assertEquals(0, env.readStart(2));
    assertEquals(2, env.readEnd(2));
  }
}
