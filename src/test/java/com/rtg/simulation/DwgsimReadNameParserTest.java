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

package com.rtg.simulation;

import junit.framework.TestCase;

/**
 */
public class DwgsimReadNameParserTest extends TestCase {

  public void testName() {
    String name = "@CFTR.3.70s_282_148_0_1_0_0_4:0:0_4:0:0_10c/1";
    String read = "CTCCTTACTGGGAAAGAATCATAGCTTCCTATGACCCCGGATAACAAGGAGGTAACGCTCTATCGCGATTTATC";
    final DwgsimReadNameParser parser = new DwgsimReadNameParser();

    assertTrue(parser.setReadInfo(name, read.length()));

    assertEquals("268", parser.readName());
    assertEquals(282, parser.templatePosition());
    assertTrue(parser.forwardFrame());
    assertEquals(0, parser.substitutions());
    assertEquals(4, parser.insertions());
    assertEquals(0, parser.deletions());
    assertEquals(4, parser.numMismatches());

    name = "@CFTR.3.70s_35_171_0_1_0_0_2:0:0_3:0:0_15b/2";
    read = "CAGAGAAATTAAATTTCATTTATTTCTCATAAAATACCCTGCAGAATTTCAAACACAAGACTTTTGCATCTTT";

    assertTrue(parser.setReadInfo(name, read.length()));
    assertEquals("347", parser.readName());
    assertEquals(171, parser.templatePosition());
    assertFalse(parser.forwardFrame());
    assertEquals(3, parser.numMismatches());

  }

}
