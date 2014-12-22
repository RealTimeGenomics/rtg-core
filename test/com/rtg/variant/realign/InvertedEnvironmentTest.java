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
public class InvertedEnvironmentTest extends TestCase {

  public void test1() {
    final Environment env = new InvertedEnvironment(new EnvironmentImplementation(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        3, //start
        new byte[] {1, 2, 3, 4, 0}, //read
        new double[]{0.1, 0.001, 0.002, 0.003, 0.25} //quality
        ));

    Exam.integrity(env);
    final String exp = ""
        + "Environment" + LS
        + "Template [7..2)" + LS
        + " [1]C" + LS
        + " [2]A" + LS
        + " [3]N" + LS
        + " [4]T" + LS
        + "+[5]G" + LS
        + "+[6]C" + LS
        + "Read" + LS
        + "[0]N 0.250" + LS
        + "[1]T 0.003" + LS
        + "[2]G 0.002" + LS
        + "[3]C 0.001" + LS
        + "[4]A 0.100" + LS
        ;
    assertEquals(exp, env.toString());
    assertEquals(5, env.readLength());
    assertEquals(7, env.absoluteTemplatePosition(0));
    assertEquals(3, env.absoluteTemplatePosition(4));

    assertEquals(0.25, env.quality(0));
    assertEquals(0.1, env.quality(4));

    assertEquals(0, env.read(0));
    assertEquals(1, env.read(4));

    assertEquals(0, env.template(-1));
    assertEquals(0, env.template(0));
    assertEquals(2, env.template(1));
    assertEquals(1, env.template(7));
    assertEquals(0, env.template(8));
  }

  public void test2() {
    final Environment env = new InvertedEnvironment(new EnvironmentImplementation(
        2, //maxShift
        new byte[] {1, 2, 3, 4, 0, 1, 2}, //template
        0, //start
        new byte[] {1, 2}, //read
        new double[]{0.1, 0.001} //quality
        ), 5);
    Exam.integrity(env);
    final String exp = ""
        + "Environment" + LS
        + "Template [6..4)" + LS
        + " [0]C" + LS
        + " [1]A" + LS
        + "+[2]N" + LS
        + "+[3]T" + LS
        + "Read" + LS
        + "[0]C 0.001" + LS
        + "[1]A 0.100" + LS
        ;
    assertEquals(exp, env.toString());
  }

  public void test3() {
    final Environment env = new InvertedEnvironment(new EnvironmentImplementation(
        0, //maxShift
        new byte[] {1}, //template
        0, //start
        new byte[] {1}, //read
        new double[]{0.1} //quality
        ));
    Exam.integrity(env);
    final String exp = ""
        + "Environment" + LS
        + "Template [0..-1)" + LS
        + " [0]A" + LS
        + "Read" + LS
        + "[0]A 0.100" + LS
        ;
    assertEquals(exp, env.toString());
  }
}
