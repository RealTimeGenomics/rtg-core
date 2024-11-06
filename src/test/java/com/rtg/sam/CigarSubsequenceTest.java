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
package com.rtg.sam;

import com.rtg.AbstractTest;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.Match;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 *
 */
public class CigarSubsequenceTest extends AbstractTest {

  public static SAMRecord makeSamRecord(final int alignStart, final String readString, final String cigar) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(alignStart);
    sam.setReadString(readString);
    sam.setCigarString(cigar);
    final byte[] qual = new byte[readString.length()];
    for (int i = 0; i < readString.length(); ++i) {
      qual[i] = (byte) (i + 1);
    }
    sam.setBaseQualities(qual);
    return sam;
  }

  //by default use CigarFormatter.cigarSubsequence to produce match. Override to test other implementations
  protected AlignmentMatch getAlignmentMatch(VariantAlignmentRecord rec, int start, int end, VariantParams params) {
    return CigarFormatter.cigarSubsequence(rec, null, start, end, params);
  }

  protected AlignmentMatch checkCigarSubSequence(final String read, final String cigar, final int start, final int end, final String expected) {
    final SAMRecord sam = makeSamRecord(2, read, cigar);
    final AlignmentMatch am = getAlignmentMatch(new VariantAlignmentRecord(sam), start + 1, end + 1, new VariantParamsBuilder().create());
    if (expected == null) {
      assertNull(am);
    } else {
      assertNotNull(am);
      assertEquals(expected, am.toString());
    }
    return am;
  }

  public void testCigarSubSequenceDeletes() {
    checkCigarSubSequence("ACGTACGTACGT", "4=4D4=", 2, 5, "GT");
    checkCigarSubSequence("ACGTACGTACGT", "4=4D4=", 4, 7, "");
    checkCigarSubSequence("ACGTACGTACGT", "4D4=4=", 2, 6, "AC");
    checkCigarSubSequence("ACGTACGTACGT", "4=4D4=", 2, 6, "GT");
  }

  public void testCigarSubSequenceSoft() {
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "12=", 0, 2, "AC"), 0, 12);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "2H1S11=", 0, 2, "CG"), 1, 12);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6=6S", 2, 6, "GTAC"), 0, 6);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6=6S", 2, 8, "GTAC~"), 0, 6);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6S6=", 1, 5, "TACG"), 6, 12);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6S6=", 0, 4, "GTAC"), 6, 12);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6S6=", -1, 4, "~GTAC"), 6, 12);

    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "6=6S", 0, 2, "AC"), 0, 6);
    checkSoftClipBounds(checkCigarSubSequence("ACGTACGTACGT", "4H12=4H", 0, 12, "ACGTACGTACGT"), 0, 12);
  }

  public void testCigarSubSequenceP() {
    checkCigarSubSequence("ACGTACGTACGT", "40P6I30P6=4P", 0, 6, "ACGTACGTACGT");
  }

  public void testFoo() {
    checkCigarSubSequence("ACGTACGTACGT", "2=4M6=", 2, 6, "GTAC");
  }


  private static void checkSoftClipBounds(AlignmentMatch match, int expectedLeft, int expectedRight) {
    assertEquals(expectedLeft, match.getSoftClipLeft());
    assertEquals(expectedRight, match.getSoftClipRight());
  }

  public void testCigarSubSequenceLeftFixed() {
    checkCigarSubSequence("ACGTACGTACGT", "4=4N4=", 2, 5, "GT~");
    checkCigarSubSequence("ACGTACGTACGT", "4=1=1N2=4=", 2, 6, "GTA~");
    checkCigarSubSequence("ACGTACGTACGT", "4=1=1N2=4=", -1, 6, "~ACGTA~");
  }

  public void testCigarSubSequenceRightFixed() {
    checkCigarSubSequence("ACGTACGTACGT", "4=1=1N2=4=", 5, 9, "~CGT");
    checkCigarSubSequence("ACGTACGTACGT", "4=3N5=", 6, 10, "~ACG");
    checkCigarSubSequence("ACGTACGTACGT", "12=", 6, 13, "GTACGT~");
  }

  public void testCigarSubSequenceNeitherFixed() {
    checkCigarSubSequence("ACGTACGTACGT", "4=1N2=1N4=", 4, 8, "~AC~");
    checkCigarSubSequence("ACGTACGTACGT", "4=2N2=1N3=", 4, 9, "~AC~");
  }

  public void testCigarSubSequenceInvalid() {
    checkCigarSubSequence("ACGTACGTACGT", "2=40I", 1, 6, null);
    checkCigarSubSequence("ACGTACGTACGT", "2=4N6I", 1, 6, null);
    checkCigarSubSequence("ACGTACGTACGT", "2=4N6=", 1, 6, "C~");
    checkCigarSubSequence("ACGTACGTACGT", "4=1N1=1N1=4=", 2, 6, null);
    checkCigarSubSequence("ACGTACGTACGT", "4=1=1N2=4=", 2, 7, null);
    checkCigarSubSequence("ACGTACGTACGT", "2=4N6=", 2, 6, null);
    checkCigarSubSequence("ACGTACGTACGT", "1=2I4N6=", 1, 7, null);
  }

  public void testCigarSubSequenceAllCigars() {
    checkCigarSubSequence("ACGTACGTACGT", "4S4=4=", 1, 4, "CGT"); //very tricky - the reference position starts at the end of the S cigars
    checkCigarSubSequence("ACGTACGTACGT", "4=4=4S2H", 6, 10, "GT~");
    checkCigarSubSequence("ACGTACGTACGT", "4=4=4=", 2, 6, "GTAC");
    checkCigarSubSequence("ACGTACGTACGT", "4=4M4=", 2, 6, "GTAC");
    checkCigarSubSequence("ACGTACGTACGT", "4=4X4=", 2, 6, "GTAC");
    checkCigarSubSequence("ACGTACGTACGT", "4=4P4=", 2, 6, "GTAC");
    checkCigarSubSequence("ACGTACGTACGT", "4=2I6=", 2, 6, "GTACGT");

    checkCigarSubSequence("ACGTACGTACGT", "4X4=4=", 2, 6, "GTAC");
    checkCigarSubSequence("ACGTACGTACGT", "4P4=4=", 2, 6, "GTAC");

    checkCigarSubSequence("ACGTACGTACGT", "2I4=4=4=", 2, 6, "ACGT");
  }

  public void testCigarSubSequenceSimple() {
    checkCigarSubSequence("ACGT", "4M", 0, 4, "ACGT");
    checkCigarSubSequence("ACGT", "4M", 0, 6, "ACGT~");
    checkCigarSubSequence("ACGT", "4M", 5, 8, null);
    checkCigarSubSequence("TTACGTAA", "2=2I4=", 2, 4, "ACGT");
    checkCigarSubSequence("TTACGTAA", "2=2I4=", 3, 5, "TA");
    checkCigarSubSequence("ACGTACGTACGT", "12M", 0, 12, "ACGTACGTACGT");
    checkCigarSubSequence("ACGT", "4X", 0, 4, "ACGT");
    checkCigarSubSequence("ACGT", "4=", 0, 4, "ACGT");
    checkCigarSubSequence("ACGT", "4M", 0, 4, "ACGT");
    checkCigarSubSequence("ACGT", "4M4S", 1, 4, "CGT");
    checkCigarSubSequence("AACGT", "1S4M", 1, 4, "CGT");
    checkCigarSubSequence("ACGT", "4M", 0, 3, "ACG");
    checkCigarSubSequence("ACGT", "2M2M", 0, 4, "ACGT");
    checkCigarSubSequence("ACGT", "2M2M", 1, 4, "CGT");
  }

  //test when counts in cigar >= 10 (take care to cover first and second loop)
  public void testCigarSubSequenceLong() {
    checkCigarSubSequence("CCACGGGGGGGGGTC", "3=10I2=", 4, 6, "C~");
    checkCigarSubSequence("CCACGGGGGGGGGTC", "3=10I2=", 3, 6, "CGGGGGGGGGTC~");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 1, 2, "C");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 2, 2, "");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 2, 3, "ACGGGGGGGGG");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 2, 4, "ACGGGGGGGGGT");
    checkCigarSubSequence("CCACGGGGGGGGGT", "2=10I2=", 2, 3, "ACGGGGGGGGG");
    checkCigarSubSequence("CCACGGGGGGGGGT", "1=10I3=", 2, 3, "G");
    checkCigarSubSequence("CCACGGGGGGGGAT", "1=10I3=", 2, 3, "A");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 3, 3, "CGGGGGGGGG");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 3, 4, "CGGGGGGGGGT");
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 4, 4, null);
    checkCigarSubSequence("CCACGGGGGGGGGT", "3=10I1=", 4, 5, null);
  }

  //test when insert at end of selected region
  public void testCigarSubSequenceInsertAtEnd() {
    checkCigarSubSequence("CCGT", "2=1I1=", 2, 2, "G");
    checkCigarSubSequence("CCGT", "3I1=", 0, 0, "CCG");
    checkCigarSubSequence("CCGT", "3S1=", 0, 0, null);

    checkCigarSubSequence("ACGTACGTACGTACGT", "8=4I4=", 0, 8, "ACGTACGTACGT");
  }

  //start and end same (and inserts there)
  public void testCigarSubSequenceSame() {
    checkCigarSubSequence("ACGT", "4=", 0, 0, null);
    checkCigarSubSequence("ACGT", "1I3=", 0, 0, "A");
  }


  //test non-null quality
  public void testCigarSubSequenceQuality1() {
    final SAMRecord sam = makeSamRecord(0, "TTACGTAA", "2=2I4=");
    final Match cs = CigarFormatter.cigarSubsequence(new VariantAlignmentRecord(sam), null, 2, 4,  new VariantParamsBuilder().create());
    assertEquals("TA", cs.toString());
    assertEquals("[0.60, 0.70, ]", cs.qualityString());
  }

  //check case when no quality scores provided
  public void testCigarSubSequenceQuality2() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(0);
    sam.setReadString("TTACGTAA");
    sam.setCigarString("2H2=2I4=3H");
//    final byte[] qual = new byte["TTACGTAA".length()];
//    for (int i = 0; i < "TTACGTAA".length(); ++i) {
//      qual[i] = (byte) (i + 1);
//    }
    //sam.setBaseQualities(qual);
    final Match cs = CigarFormatter.cigarSubsequence(new VariantAlignmentRecord(sam), null, 2, 3,  new VariantParamsBuilder().create());
    assertEquals("T", cs.toString());
    assertEquals("[2.00, ]", cs.qualityString());
  }

}
