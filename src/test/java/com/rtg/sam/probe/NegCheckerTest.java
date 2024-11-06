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

package com.rtg.sam.probe;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 *
 */
public class NegCheckerTest extends TestCase {

  public void test() {
    final NegChecker neg = new NegChecker(10);
    final SAMRecord rec = PosCheckerTest.createRecord("AGGTTTGG", "1=1X1=1X3=1X");
    neg.setAlignmentEnd(rec, null, 1006);
    assertEquals("AGGTTT", rec.getReadString());
    assertEquals("1=1X1=1X2=", rec.getCigarString());
  }

  public void testTrimAllProbe() {
    final NegChecker neg = new NegChecker(10);
    final SAMRecord rec = PosCheckerTest.createRecord("AGGTTTGG", "1=1X1=1X3=1X");
    neg.setAlignmentEnd(rec, null, 990);
    assertEquals("*", rec.getReadString());
    assertEquals("*", rec.getCigarString());
  }
}
