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
package com.rtg.mode;

import java.util.Arrays;
import java.util.Locale;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class DNATest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.DNA()}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(DNA.class, "[N, A, C, G, T]");
    try {
      DNA.valueOf('u');
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("u", e.getMessage());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.DNA#ignore()}.
   */
  public final void testIsUnknown() {
    assertTrue(DNA.N.ignore());
    assertEquals(0, DNA.N.ordinal());
    assertTrue(DNA.valueOf("N").ignore());
    for (final DNA dna : DNA.values()) {
      if (dna != DNA.N) {
        assertFalse(dna.ignore());
      }
      if (dna.ignore()) {
        assertTrue(dna.ordinal() < SequenceType.DNA.firstValid());
      } else {
        assertTrue(dna.ordinal() >= SequenceType.DNA.firstValid());
      }
    }
  }

  /**
   * Test method for {@link com.rtg.mode.DNA#type()}.
   */
  public final void testType() {
    for (final DNA dna : DNA.values()) {
      assertEquals(SequenceType.DNA, dna.type());
    }
  }

  public final void testChars() {
    final char[] chars = DNA.valueChars();
    assertEquals('N', chars[0]);
    assertEquals('A', chars[1]);
    assertEquals('C', chars[2]);
    assertEquals('G', chars[3]);
    assertEquals('T', chars[4]);
    assertEquals(5, chars.length);
  }

  /**
   * Test method for {@link com.rtg.mode.DNA#complement()}.
   */
  public final void testComplement() {
    for (final DNA dna : DNA.values()) {
      assertEquals(dna, dna.complement().complement());
    }
    assertEquals(DNA.N, DNA.N.complement());
    assertEquals(DNA.T, DNA.A.complement());
    assertEquals(DNA.G, DNA.C.complement());
    assertEquals(DNA.C, DNA.G.complement());
    assertEquals(DNA.A, DNA.T.complement());
    assertEquals(DNA.N.ordinal(), DNA.complement((byte) DNA.N.ordinal()));
    assertEquals(DNA.A.ordinal(), DNA.complement((byte) DNA.T.ordinal()));
  }

  /**
   * Test method for {@link com.rtg.mode.DNA}.
   */
  public final void testCodeComplement() {
    for (final DNA dna : DNA.values()) {
      assertEquals(dna.complement().ordinal(), DNA.COMPLEMENT[dna.ordinal()]);
    }
  }


  /**
   * Test method for {@link com.rtg.mode.DNA#toString()}.
   */
  public final void testToString() {
    for (final DNA dna : DNA.values()) {
      final String str = dna.toString();
      assertEquals(1, str.length());
      assertEquals(str.toUpperCase(Locale.getDefault()), str);
    }
  }

  public void testName() {
    assertEquals("T", DNA.T.name());
    assertEquals("N", DNA.N.name());
  }

  public void testComplementInPlace() {
    final byte[] nt = DNA.stringDNAtoByte("NACGTNACGT");
    assertEquals("[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]", Arrays.toString(nt));
    DNA.reverseComplementInPlace(nt, 2, 5);
    assertEquals("[0, 1, 4, 0, 1, 2, 3, 2, 3, 4]", Arrays.toString(nt));
    DNA.complementInPlace(nt, 1, 2);
    assertEquals("[0, 4, 1, 0, 1, 2, 3, 2, 3, 4]", Arrays.toString(nt));

    final byte[] nt2 = DNA.stringDNAtoByte("ACGTNTCGA");
    assertEquals("[1, 2, 3, 4, 0, 4, 2, 3, 1]", Arrays.toString(nt2));
    DNA.reverseComplementInPlace(nt2);
    assertEquals("[4, 2, 3, 1, 0, 1, 2, 3, 4]", Arrays.toString(nt2));
  }

  public void testGetDNA() {
    assertEquals(0, DNA.getDNA('N'));
    assertEquals(0, DNA.getDNA('n'));
    assertEquals(1, DNA.getDNA('A'));
    assertEquals(1, DNA.getDNA('a'));
    assertEquals(2, DNA.getDNA('C'));
    assertEquals(2, DNA.getDNA('c'));
    assertEquals(3, DNA.getDNA('G'));
    assertEquals(3, DNA.getDNA('g'));
    assertEquals(4, DNA.getDNA('T'));
    assertEquals(4, DNA.getDNA('t'));
  }

  public void testPopulateComplementArray() {
    final byte[] comps = new byte[1];
    final DNA[] dna = {new DNA(0, DNASimple.N) {
      @Override
      public DNA complement() {
        return DNA.A;
      }
    }};
    DNA.populateComplementArray(comps, dna);
    assertEquals(DNA.A.ordinal(), comps[0]);
  }

  public void testStringDNAtoBytes() {
    assertEquals("[0, 1, 2, 3, 4]", Arrays.toString(DNA.stringDNAtoByte("NACGT")));
    assertEquals("[1, 1, 2, 3, 4]", Arrays.toString(DNA.stringDNAtoByte("AACGT")));
    assertEquals("[]", Arrays.toString(DNA.stringDNAtoByte("")));
  }

  public void testByteDNAtoByte() {
    assertEquals("[0, 1, 2, 3, 4]", Arrays.toString(DNA.byteDNAtoByte("NACGT".getBytes())));
    assertEquals("[1, 1, 2, 3, 4]", Arrays.toString(DNA.byteDNAtoByte("AACGT".getBytes())));
    assertEquals("[]", Arrays.toString(DNA.byteDNAtoByte("".getBytes())));
  }
}

