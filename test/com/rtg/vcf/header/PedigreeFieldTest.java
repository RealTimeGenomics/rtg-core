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
public class PedigreeFieldTest extends TestCase {

  public void testSomeMethod() {
    final String io = "##PEDIGREE=<Child=NA19240,Mother=NA19239,Father=NA19238>";
    final PedigreeField f = new PedigreeField(io);
    assertEquals("NA19240", f.getChild());
    assertEquals("NA19239", f.getMother());
    assertEquals("NA19238", f.getFather());
    assertEquals(io, f.toString());
  }
}
