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
package com.rtg.vcf;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ChildPhasingVcfAnnotatorTest extends TestCase {


  public void testPhaseCall() {
    // Mendelian calls
    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/0", "."));
    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("1/0", ".", "1/0"));
    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "0/1", "0/1"));

    assertEquals("0|0", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/0", "0/0"));
    assertEquals("0|0", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", ".", "0/0"));
    assertEquals("0|0", ChildPhasingVcfAnnotator.phaseDiploidCall(".", ".", "0/0"));

    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "1", "2/1"));
    assertEquals("0|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "1", "0/1"));

    assertEquals("0|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/1", "0/1"));
    assertEquals("0|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/1", "1/0"));
    assertEquals("1|0", ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "0/0", "0/1"));
    assertEquals("1|0", ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "0/2", "1/0"));
    assertEquals("1|0", ChildPhasingVcfAnnotator.phaseDiploidCall("1/1", "0/2", "0/1"));
    assertEquals("1|2", ChildPhasingVcfAnnotator.phaseDiploidCall("1/1", "0/2", "1/2"));
    assertEquals("1|2", ChildPhasingVcfAnnotator.phaseDiploidCall("1/1", "2/2", "1/2"));
    assertEquals("1|0", ChildPhasingVcfAnnotator.phaseDiploidCall("1/1", ".", "1/0"));

    // Non-mendelian calls
    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/0", "0/1"));
    assertNull(ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "1", "2/1"));
    assertEquals("2|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "1", "2/1"));
    assertEquals("1|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "0/0", "1/1"));
    assertEquals("2|1", ChildPhasingVcfAnnotator.phaseDiploidCall("0/0", "1/0", "1/2"));
    assertEquals("1|2", ChildPhasingVcfAnnotator.phaseDiploidCall("0/1", "0/0", "1/2"));
  }

}
