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
package com.rtg.variant;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
// TODO dramatically improve the quality of the testing
public class SamToMatchCigarTest extends TestCase {

  private static class MockCigarParser implements ReadParserInterface {
    private final StringBuilder mSb = new StringBuilder();

    @Override
    public void toMatcher(AbstractMachineErrorParams me, VariantAlignmentRecord var, int qdefault, final byte[] templateBytes) {
      mSb.append("called");
    }

    @Override
    public String toString() {
      return mSb.toString();
    }
  }

  public void test() {
    // all this really does is set up the contexts and prove that the call is
    // made - needs to be made more thorough
    final VariantParams params = VariantParams.builder().create();
    final ReadParserInterface parser = new MockCigarParser();
    final SamToMatch stm = new SamToMatchCigar(params, parser, new DefaultMachineErrorChooser());
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadBases("ACGTN".getBytes());
    sam.setCigarString("4=");
    sam.setAlignmentStart(42);
    final VariantAlignmentRecord rec = new VariantAlignmentRecord(sam);
    assertTrue(stm.process(new byte[] {0, 0}, rec));
    assertEquals("called", parser.toString());
    assertEquals(41, stm.start(rec));
  }
}
