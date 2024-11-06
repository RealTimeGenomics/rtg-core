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
package com.rtg.variant;

import com.rtg.util.machine.MachineType;

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
    public void toMatcher(VariantAlignmentRecord var, MachineType machineType, int qdefault, final byte[] templateBytes) {
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
