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

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class EnvironmentImplementationTest extends TestCase {

  protected Environment createEnv(int maxShift, byte[] template, int start, byte[] read, double[] quality) {
    return new EnvironmentImplementation(maxShift, template, start, read, quality);
  }

  public void test1() {
    final Environment env = createEnv(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        3, //start
        new byte[] {1, 2, 3, 4, 0}, //read
        new double[]{0.1, 0.001, 0.002, 0.003, 0.25} //quality
        );
    Exam.integrity(env);
    Exam.globalIntegrity(env);
    assertEquals(2, env.maxShift());
    final String exp = ""
        + "Environment" + LS
        + "Template [3..8)" + LS
        + "-[-2]C" + LS
        + "-[-1]G" + LS
        + " [0]T" + LS
        + " [1]N" + LS
        + " [2]A" + LS
        + " [3]C" + LS
        + "Read" + LS
        + "[0]A 0.100" + LS
        + "[1]C 0.001" + LS
        + "[2]G 0.002" + LS
        + "[3]T 0.003" + LS
        + "[4]N 0.250" + LS
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
  }

  public void test2() {
    final Environment env = createEnv(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        0, //start
        new byte[] {1, 2}, //read
        new double[]{0.1, 0.001} //quality
        );
    Exam.integrity(env);
    final String exp = ""
        + "Environment" + LS
        + "Template [0..2)" + LS
        + " [0]A" + LS
        + " [1]C" + LS
        + "+[2]G" + LS
        + "+[3]T" + LS
        + "Read" + LS
        + "[0]A 0.100" + LS
        + "[1]C 0.001" + LS
        ;
    assertEquals(exp, env.toString());
  }

  public void test3() {
    final Environment env = createEnv(
        0, //maxShift
        new byte[] {1}, //template
        0, //start
        new byte[] {1}, //read
        new double[]{0.1} //quality
        );
    Exam.integrity(env);
    final String exp = ""
        + "Environment" + LS
        + "Template [0..1)" + LS
        + " [0]A" + LS
        + "Read" + LS
        + "[0]A 0.100" + LS
        ;
    assertEquals(0, env.maxShift());
    assertEquals(1, env.templateLength());
    assertEquals(exp, env.toString());
  }

  public void test4() {
    final Environment env = createEnv(
        0, //maxShift
        new byte[] {1}, //template
        0, //start
        new byte[] {1}, //read
        null //quality
        );
    Exam.integrity(env);
    assertEquals(0.01, env.quality(1));
  }
}
