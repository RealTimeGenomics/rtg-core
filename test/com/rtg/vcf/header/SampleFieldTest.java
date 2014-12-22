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
package com.rtg.vcf.header;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SampleFieldTest extends TestCase {

  public void testSomeMethod() {
    final String io = "##SAMPLE=<ID=yo, Genomes=cancer;normal, Mixture=roobar, Description=\"fun for the whole family\">";
    final SampleField f = new SampleField(io);
    assertEquals(io.replaceAll(",\\s+", ","), f.toString());
    assertEquals("yo", f.getId());
    assertEquals("cancer;normal", f.getGenomes());
    assertEquals("roobar", f.getMixture());
    assertEquals("fun for the whole family", f.getDescription());

    final String io2 = "##SAMPLE=<ID=%$&*,Genomes=A;B,Mixture=robbyy,Description=\"oh rearry\">";
    final SampleField f2 = new SampleField(io2);
    assertEquals(io2, f2.toString());
    assertEquals("%$&*", f2.getId());
    assertEquals("A;B", f2.getGenomes());
    assertEquals("robbyy", f2.getMixture());
    assertEquals("oh rearry", f2.getDescription());
  }
}
