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

import java.util.Locale;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class SuperCigarParserTest extends TestCase {

  static class SuperCigarParserTester extends SuperCigarParser {

    StringBuilder mChunks = new StringBuilder();
    StringBuilder mTemplateString = new StringBuilder();
    StringBuilder mReadString = new StringBuilder();
    private static final char[] DNA_CHAR = DNA.valueChars();
    byte[] mReadSeq = null;

    @Override
    protected void doChunk(char cmd, int size) throws BadSuperCigarException {
      super.doChunk(cmd, size);
      mChunks.append(size);
      mChunks.append(cmd);
      mChunks.append(" ");
    }

    @Override
    protected void doReadHardClip() {
      mTemplateString.append(DNA.N);
    }

    @Override
    protected void doReadSoftClip(int readNt) {
      mReadString.append(DNA_CHAR[readNt]);
      mTemplateString.append(DNA.N);
      Exam.assertEquals("read[" + getReadPosition() + "]", readNt, mReadSeq[getReadPosition()]);
    }

    @Override
    protected void doTemplateOverlap() {
      mTemplateString.setLength(mTemplateString.length() - 1);
    }

    @Override
    protected void doTemplateSkip(int templateNt) {
      mTemplateString.append('.');
      Exam.assertEquals(templateNt, templateNt(getTemplatePosition()));
    }

    @Override
    protected void doReadOnly(int readNt) {
      mReadString.append(DNA_CHAR[readNt]);
      Exam.assertEquals(readNt, mReadSeq[getReadPosition()]);
    }

    @Override
    protected void doTemplateOnly(int templateNt) {
      mTemplateString.append(DNA_CHAR[templateNt]);
      Exam.assertEquals(templateNt, templateNt(getTemplatePosition()));
    }

    @Override
    protected void doSubstitution(int readNt, int templateNt) {
      mReadString.append(DNA_CHAR[readNt]);
      mTemplateString.append(DNA_CHAR[templateNt]);
      Exam.assertEquals(templateNt, templateNt(getTemplatePosition()));
      Exam.assertEquals(readNt, mReadSeq[getReadPosition()]);
    }

    @Override
    protected void doEquality(int readNt, int nt) {
      mReadString.append(DNA_CHAR[nt]);
      mTemplateString.append(DNA_CHAR[nt]);
      Exam.assertEquals(nt, templateNt(getTemplatePosition()));
      Exam.assertEquals(nt, mReadSeq[getReadPosition()]);
    }
  }

  public void testParse() throws BadSuperCigarException {
    final String template = " ccTTT    cccccTTTTTTTTTTTTTTTTTTTT     ttt".toUpperCase(Locale.ROOT).replaceAll(" ", "");
    final String read =     "G  CCCGGGG     TTTTTTTTTTTTTTTTTTTT A TTTTT".replaceAll(" ", "");
    final SuperCigarParserTester parser = new SuperCigarParserTester();
    parser.mReadSeq = DnaUtils.encodeString(read);
    parser.setTemplate(DnaUtils.encodeString(template));
    parser.setTemplateStart(0);
    assertEquals(0, parser.getTemplateStartPosition());
    parser.setCigar("3H1S2D3X4I5N20=2B1I5=2H", "gCCCggggA");
    parser.parse();
    assertEquals("3H 1S 2D 3X 4I 5N 20= 2B 1I 5= 2H ", parser.mChunks.toString());
    assertEquals(read, parser.mReadString.toString());
    assertEquals("NNNN" + template + "NN", parser.mTemplateString.toString().replaceAll("\\.", "C"));
    //assertEquals(-1, parser.getReadOverlapStart());
    /*
    // now do the same thing in reverse
    parser.mChunks = new StringBuilder();
    parser.mReadString = new StringBuilder();
    parser.mTemplateString = new StringBuilder();
    parser.parseReverse();
    assertEquals("2H 5= 1I 2B 20= 5N 4I 3X 2D 1S 3H ", parser.mChunks.toString());
    assertEquals(read, parser.mReadString.reverse().toString());
    assertEquals("NNN" + template + "NN", parser.mTemplateString.reverse().toString().replaceAll("\\.", "C"));
    */
  }

  public void testStandardCigar() throws BadSuperCigarException {
    final String template = " ccTTT    cccccTTTTTTTTTTTTTTTTTT   ggg".toUpperCase(Locale.ROOT).replaceAll(" ", "");
    final String read =     "G  CCCGGGG     TTTTTTTTTTTTTTTTTT A GGG".replaceAll(" ", "");
    final SuperCigarParserTester parser = new SuperCigarParserTester();
    parser.mReadSeq = DnaUtils.encodeString(read);
    parser.setTemplate(DnaUtils.encodeString(template));
    parser.setTemplateStart(0);
    assertEquals(0, parser.getTemplateStartPosition());
    parser.setStandardCigar("3H1S2D3X4I5N18=1I3=2H", read);
    parser.parse();
    assertEquals("3H 1S 2D 3X 4I 5N 18= 1I 3= 2H ", parser.mChunks.toString());
    assertEquals(read, parser.mReadString.toString());
    assertEquals("NNNN" + template + "NN", parser.mTemplateString.toString().replaceAll("\\.", "C"));
  }

  public void testSuperWipesRead() {
    final SuperCigarParserTester parser = new SuperCigarParserTester();
    parser.setStandardCigar("3H1S2D3X4I5N18=1I3=2H", new byte[] {0, 3, 2}, 3);
    assertEquals(3, parser.getReadLength());
    parser.setCigar("3H", "NNAN");
    assertEquals(0, parser.getReadLength());
  }

  public void testIllegalStart() {
    final SuperCigarParser parser = new SuperCigarParser();
    parser.setCigar("4M", "TTTT");
    parser.setTemplate(new byte[] {1, 3});
    parser.setTemplateStart(0);
    try {
      parser.parse();
      fail();
    } catch (final BadSuperCigarException e) {
      assertEquals("bad cigar start=0, pos=2, len=2, cigar 4M", e.getMessage());
    }
  }

  public void testIllegalCigar() {
    final SuperCigarParser parser = new SuperCigarParser();
    parser.setCigar("1P", "CCC");
    parser.setTemplate(new byte[] {1, 3});
    parser.setTemplateStart(0);
    try {
      parser.parse();
      fail();
    } catch (final BadSuperCigarException e) {
      assertEquals("Illegal command in supercigar: 1P", e.getMessage());
    }
  }

  /**
   * Test overlap size calculations.  These are independent of left/right arm.
   * @throws BadSuperCigarException on error.
   */
  public void testOverlap() throws BadSuperCigarException {
    final SuperCigarParser parser = new SuperCigarParser();
    // the simplest case
    parser.setCigar("30=1B5=", null);
    parser.setTemplate(new byte[135]); // all N
    parser.setTemplateStart(100);
    assertEquals(-1, parser.getReadOverlapStart());
    assertEquals(-1, parser.getReadOverlapMiddle());
    assertEquals(-1, parser.getReadOverlapEnd());
    assertEquals(-1, parser.getTemplateOverlapStart());
    assertEquals(-1, parser.getTemplateOverlapEnd());
    parser.parse();
    assertEquals(29, parser.getReadOverlapStart());
    assertEquals(30, parser.getReadOverlapMiddle());
    assertEquals(31, parser.getReadOverlapEnd());
    assertEquals(129, parser.getTemplateOverlapStart());
    assertEquals(130, parser.getTemplateOverlapEnd());

    // a cigar with complex stuff in the second half of the overlap
    parser.setCigar("5=3B1=1X1I1=2I25=", "CGTT");
    parser.setTemplateStart(100);
    assertEquals(-1, parser.getReadOverlapStart());
    assertEquals(-1, parser.getReadOverlapMiddle());
    assertEquals(-1, parser.getReadOverlapEnd());
    assertEquals(-1, parser.getTemplateOverlapStart());
    assertEquals(-1, parser.getTemplateOverlapEnd());
    parser.parse();
    assertEquals(2, parser.getReadOverlapStart());
    assertEquals(5, parser.getReadOverlapMiddle());
    assertEquals(9, parser.getReadOverlapEnd());
    assertEquals(102, parser.getTemplateOverlapStart());
    assertEquals(105, parser.getTemplateOverlapEnd());

    // and another cigar with complex stuff in the left half of the overlap
    parser.setCigar("1S1=1I1D1X1=3B1I2=2D15=8N10=", "ACGT");
    parser.setTemplateStart(0);
    assertEquals(-1, parser.getReadOverlapStart());
    assertEquals(-1, parser.getReadOverlapMiddle());
    assertEquals(-1, parser.getReadOverlapEnd());
    assertEquals(-1, parser.getTemplateOverlapStart());
    assertEquals(-1, parser.getTemplateOverlapEnd());
    parser.parse();
    assertEquals(3, parser.getReadOverlapStart());
    assertEquals(5, parser.getReadOverlapMiddle());
    assertEquals(8, parser.getReadOverlapEnd());
    assertEquals(1, parser.getTemplateOverlapStart());
    assertEquals(4, parser.getTemplateOverlapEnd());

    // This one the left side of the overlap goes right back to 0.
    // In one run, this returned a stale value and caused an array bounds error.
    parser.setCigar("1=1I3=4B20=7N10=", "A");
    parser.setTemplateStart(0);
    assertEquals(-1, parser.getReadOverlapStart());
    assertEquals(-1, parser.getReadOverlapMiddle());
    assertEquals(-1, parser.getReadOverlapEnd());
    assertEquals(-1, parser.getTemplateOverlapStart());
    assertEquals(-1, parser.getTemplateOverlapEnd());
    parser.parse();
    assertEquals(0, parser.getReadOverlapStart());
    assertEquals(5, parser.getReadOverlapMiddle());
    assertEquals(9, parser.getReadOverlapEnd());
    assertEquals(0, parser.getTemplateOverlapStart());
    assertEquals(4, parser.getTemplateOverlapEnd());

    /*    0    5    0    5    0    5    0    5      0    5
    final TTCGAGTATGGAGCCAGTCCTCATAAGTCGGCCCCT  TCCTATGTCC
            CGAGTATGGA      CCTCATAAGTCGGCCCCTGGT
            0    5          0    5    0    5    0
            CGAGTATGGA      CCTCATAAGTCGGCCCCTTCCTGGT
                                            ****  II
            */
    parser.setCigar("10=6N20=4B2=2I1=", "GG");
    parser.setTemplateStart(2);
    assertEquals(-1, parser.getReadOverlapStart());
    assertEquals(-1, parser.getReadOverlapMiddle());
    assertEquals(-1, parser.getReadOverlapEnd());
    assertEquals(-1, parser.getTemplateOverlapStart());
    assertEquals(-1, parser.getTemplateOverlapEnd());
    parser.parse();
    assertEquals(26, parser.getReadOverlapStart());
    assertEquals(30, parser.getReadOverlapMiddle());
    assertEquals(35, parser.getReadOverlapEnd());
    assertEquals(34, parser.getTemplateOverlapStart());
    assertEquals(38, parser.getTemplateOverlapEnd());
  }

  public void testUnknowns() throws Exception {

    final SuperCigarParser parser = new SuperCigarParser() {
      @Override
      protected void doUnknownOnTemplate(int readNt, int templateNt) {
        assertEquals(DNA.getDNA('G'), readNt);
        assertEquals(12, mReadPos);
        assertEquals(1, templateNt);
      }

      @Override
      protected void doUnknownOnRead() {
        assertEquals(10, mReadPos);
      }

      @Override
      protected void doSubstitution(int readNt, int templateNt) {
        assertEquals(DNA.getDNA('T'), readNt);
        assertEquals(DNA.getDNA('A'), templateNt);
        assertEquals(11, mReadPos);
      }
    };
    // the simplest case
    parser.setCigar("10=1R1X1T5=", "TG");
    parser.setTemplate(new byte[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}); // all A
    parser.setTemplateStart(0);

    parser.parse();
  }
}
