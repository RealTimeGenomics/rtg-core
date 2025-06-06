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
package com.rtg.visualization;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class SamAssistanceSimpleTest extends TestCase {

  public final void testCigarToReadsWithInserts() throws BadSuperCigarException {
    check("ACGT", "ACGT", "4M", new int[] {0, 0, 0, 0});
    check(" ACGT", "ACGT", "4M", new int[] {1, 0, 0, 0}, 1, true);

    try {
      check("ACGT", "ACGT", "4Z", new int[] {0, 0, 0, 0});
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("Invalid cigar: 4Z", e.getMessage());
    }

    check("A  C", "ACGT", "1M2N1M", new int[] {0, 0, 0, 0});
    check(" A _ C", "ACGT", "1M2N1M", new int[] {1, 0, 1, 0}, 1, true);
    check("   A__  C", "ACGT", "1M2N1M", new int[] {2, 0, 2, 0, 0, 0, 0}, 3, true);

    check("A--C", "ACGT", "1M2D1M", new int[] {0, 0, 0, 0});
    check("....", "ACGT", "4=", "ACGT", 0, true, false);
    check("ACGT", "ACGT", "4=", "NNNN", 0, true, false);
    check("..............", "ACGTACGTACGTAC", "14=", "ACGTACGTACGTAC", 0, true, false);
    check("..N.", "ACNT", "4=", "ACGT", 0, true, false);
    check("ACGT", "ACGT", "4X", new int[] {0, 0, 0, 0});

    check("....ACGT", "ACGTACGT", "4=4I", "ACGTA", 0, true, false);

    //test "S"
    check("CG.", "ACGT", "1S2X1=", "TTT", 0, true, false);
    check("CG.", "ACGT", "1S2X1=", "TTT", 0, true, true);
    check("aCG.", "ACGT", "1S2X1=", "TTTT", 1, true, true);
    check("ACNT", "ACNT", "4X", new int[] {0, 0, 0, 0});

    check("AC_GT", "ACGT", "4M", new int[] {0, 0, 1, 0});
    check(".._..", "ACGT", "4=", "AC_GT", 0, true, false);
    check("AC_GT", "ACGT", "4X", new int[] {0, 0, 1, 0});
    check("ACG____T", "ACGT", "2M1I1M", new int[] {0, 0, 5, 0});

    check("..GTGT", "ACGTGT", "2=1I1X1I1X", "AC_A_A", 0, true, false);

    check(".-CG", "ACGT", "1=1D2X", "ACTT", 0, true, false);

    check("  .-CG", "ACGT", "1=1D2X", "TTACTT", 2, true, false);
    check("tt.-CG", "TTTACGT", "3S1=1D2X", "TTACTT", 2, true, true);
    check("tt.-CG", "TTACGT", "2S1=1D2X", "TTACTT", 2, true, true);
    check(" t.-CG", "TACGT", "1S1=1D2X", "TTACTT", 2, true, true);
    check(" t.-CGtttt", "TACGTTTTT", "1S1=1D2X4S", "TTACTTA", 2, true, true);
    check(" t.-CGttt", "TACGTTTT", "1S1=1D2X3S", "TTACTTA", 2, true, true);
    check(" t.-CGtt", "TACGTTTT", "1S1=1D2X2S", "TTACTTA", 2, true, true);
    check(" t.-CGt", "TACGTTTT", "1S1=1D2X1S", "TTACTTA", 2, true, true);

    //insert before start of read
    check("   .-CGT", "ACGT", "1=1D3X", "___AAAAA", 3, true, false);

    check("A-CG", "ACGT", "1=1D2X", new int[] {0, 0, 0, 0}, 0, false);
  }

  public void testInsertAtStart() {
    final String cigar = "2I3M";
    final String read = "TCAGT";
    final String template = "GGAGTGAAAGGGGAAAAGGGGAAAGAGAGAGA";
    final int readStart = 2;

    final SamAssistanceSimple sas = new SamAssistanceSimple();

    String s = sas.cigarToReads(cigar, read, template, readStart, true, false);
    assertEquals("TC...", s);

    s = sas.cigarToReads("2I5M", "TTGTGAA", template, 3, true, false);
    assertEquals(" TT.....", s);
  }

  public static SAMRecord makeSamRecord(final String readString, final String cigar) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadString(readString);
    sam.setCigarString(cigar);
    sam.setBaseQualities(new byte[readString.length()]);
    return sam;
  }

  void check(final String expected, final String read, final String cigar, final int[] inserts) throws BadSuperCigarException {
    check(expected, read, cigar, inserts, 0, true);
  }

  void check(final String expected, final String read, final String cigar, final int[] inserts, final int readStart, final boolean displayDots) throws BadSuperCigarException {
    final String template = insertString(inserts);
    check(expected, read, cigar, template, readStart, displayDots, false);
  }

  static String insertString(final int[] inserts) {
    final StringBuilder sb = new StringBuilder();
    for (int insert : inserts) {
      for (int j = 0; j < insert; ++j) {
        sb.append("_");
      }
      sb.append("X");
    }
    return sb.toString();
  }

  void check(final String expected, final String read, final String cigar, final String template, final int readStart, final boolean displayDots, boolean displaySoftClip) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceSimple();
    final SAMRecord sam = makeSamRecord(read, cigar);
    sam.setAlignmentStart(1);
    //System.err.println(template);
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), readStart, displayDots, displaySoftClip);
    assertEquals(1, actual.length);
    assertEquals(expected, actual[0]);
  }
}
