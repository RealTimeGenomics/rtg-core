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


/**
 */
public class EnvironmentSubreadTest extends EnvironmentImplementationTest {

  @Override
  protected Environment createEnv(int maxShift, byte[] template, int start, byte[] read, double[] quality) {
    return new EnvironmentSubread(maxShift, template, start, read, quality);
  }

  protected Environment createEnv(int maxShift, byte[] template, int start, byte[] read, int readStart, int readEnd, double[] quality) {
    return new EnvironmentSubread(maxShift, template, start, read, readStart, readEnd, quality);
  }

  public void test1Sub() {
    final Environment env = createEnv(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        3, //start
        new byte[] {1, 2, 3, 4, 0}, //read
        2,
        4,
        new double[]{0.1, 0.001, 0.002, 0.003, 0.25} //quality
        );
    Exam.integrity(env);
    Exam.globalIntegrity(env);
    assertEquals(2, env.maxShift());
    final String exp = ""
        + "Environment" + LS
        + "Template [3..5)" + LS
        + "-[-2]C" + LS
        + "-[-1]G" + LS
        + " [0]T" + LS
        + " [1]N" + LS
        + "+[2]A" + LS
        + "+[3]C" + LS
        + "Read" + LS
        + "[0]G 0.002" + LS
        + "[1]T 0.003" + LS
        ;
    assertEquals(exp, env.toString());
    assertEquals(2, env.readLength());
    assertEquals(3, env.absoluteTemplatePosition(0));
    assertEquals(7, env.absoluteTemplatePosition(4));

    assertEquals(0.002, env.quality(0));
    assertEquals(0.003, env.quality(1));

    assertEquals(3, env.read(0));
    assertEquals(4, env.read(1));

    assertEquals(7, env.templateLength());
    assertEquals(0, env.template(-4));
    assertEquals(1, env.template(-3));
    assertEquals(4, env.template(0));
    assertEquals(2, env.template(3));
    assertEquals(0, env.template(4));
  }

  //  @Override
  //  public void test2() {
  //    final Environment env = createEnv(
  //        2, //maxShift
  //        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
  //        0, //start
  //        new byte[] {1, 2}, //read
  //        new double[]{0.1, 0.001} //quality
  //        );
  //    Assert.integrity(env);
  //    final String exp = ""
  //        + "Environment" + LS
  //        + "Template 0..1" + LS
  //        + "[0]A" + LS
  //        + "[1]C" + LS
  //        + "[2]G" + LS
  //        + "[3]T" + LS
  //        + "Read" + LS
  //        + "[0]A 0.100" + LS
  //        + "[1]C 0.001" + LS
  //        ;
  //    assertEquals(exp, env.toString());
  //  }
  //
  //  @Override
  //  public void test3() {
  //    final Environment env = createEnv(
  //        0, //maxShift
  //        new byte[] {1}, //template
  //        0, //start
  //        new byte[] {1}, //read
  //        new double[]{0.1} //quality
  //        );
  //    Assert.integrity(env);
  //    final String exp = ""
  //        + "Environment" + LS
  //        + "Template 0..0" + LS
  //        + "[0]A" + LS
  //        + "Read" + LS
  //        + "[0]A 0.100" + LS
  //        ;
  //    assertEquals(0, env.maxShift());
  //    assertEquals(1, env.templateLength());
  //    assertEquals(exp, env.toString());
  //  }
  //
  //  @Override
  //  public void test4() {
  //    final Environment env = createEnv(
  //        0, //maxShift
  //        new byte[] {1}, //template
  //        0, //start
  //        new byte[] {1}, //read
  //        null //quality
  //        );
  //    Assert.integrity(env);
  //    assertEquals(0.01, env.quality(1));
  //  }
}
