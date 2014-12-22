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
