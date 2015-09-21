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
package com.rtg.index.hash.ngs.protein;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.StringWriter;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallMock;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallMock;
import com.rtg.index.hash.ngs.instances.SubstituteIndel;
import com.rtg.index.hash.ngs.instances.SubstituteProtein;
import com.rtg.util.test.RandomDna;
import com.rtg.mode.Protein;
import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;


/**
 */
public class ProteinMaskTest extends TestCase {

  protected void check(final Skeleton sk, final String protein, final String expected) throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final HashFunctionFactory f = ProteinMask.factory(sk);
    final NgsHashFunction hf = f.create(rcall, call);
    //System.err.println("NgsHashFunction:" + hf);
    Exam.integrity(hf);
    sb.append(hf.toString()).append(LS).append("number windows=").append(String.valueOf(hf.numberWindows())).append(LS);
    hf.templateSet(1234, 23);
    encode(hf, protein);
    hf.readAll(5, false);
    hf.templateForward(7);
    assertEquals(expected, sb.toString());
  }

  protected void check(final Skeleton sk) throws IOException {
    // We could use random protein here, but DNA is a subset of Protein...
    final String protein = RandomDna.random(sk.readLength());
    final HashFunctionFactory f = ProteinMask.factory(sk);
    final int indel = sk.indels();
    final SubstituteProtein sub = new SubstituteProtein(protein, f, true);
    sub.substituteProtected(sk.substitutions());
    if (indel > 0) {
      SubstituteIndel.checkIndel(f, protein, indel, 0/*cg*/);
    }
  }

  public static void encode(final NgsHashFunction hf, final String protein) {
    final Protein[] aa = Protein.values();
    for (int i = 0; i < protein.length(); i++) {
      boolean found = false;
      for (Protein anAa : aa) {
        if (anAa.name().charAt(0) == protein.charAt(i)) {
          hf.hashStep((byte) anAa.ordinal());
          found = true;
        }
      }
      if (!found) {
        throw new RuntimeException("Unknown amino acid: " + protein.charAt(i));
      }
    }
  }


  private static final String SUFFIX = "done" + LS;

  private static final String EXPECTED_1 = ""
    + "Mask l=4 w=4 s=0 i=0" + LS
    + "number windows=1" + LS
    + "set name=1234 length=23" + LS
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000001:00100100:10001000" + LS
    + "templateCall position=7 index=0" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000001:00100100:10001000" + LS
    + SUFFIX
    ;

  public void test1() throws IOException {
    final Skeleton sk = new Skeleton(4, 4, 0, 0, 1);
    check(sk);
    final String dna = "RNEP";  // the ordinals are 3(=1+2),4,8,16
    assertEquals(4, dna.length());
    check(sk, dna, EXPECTED_1);
  }

  public void test() {
    final ProteinMask hf = dummyMask();
    assertFalse(hf.isValid());
    hf.integrity(0, 0L, 0L, 0L, 0L, 0L);
    hf.hashStep((byte) 0);
    hf.integrity(1, 0L, 0L, 0L, 0L, 0L);
    hf.hashStep((byte) 1);
    hf.integrity(2, 0L, 0L, 0L, 0L, 1L);
    hf.hashStep((byte) 2);
    hf.integrity(3, 0L, 0L, 0L, 1L, 2L);
    hf.hashStep((byte) 3);
    hf.integrity(4, 0L, 0L, 0L, 3L, 5L);
    hf.hashStep((byte) 21);
    hf.integrity(5, 1L, 0L, 1L, 6L, 11L);
    hf.hashStep((byte) 8);
    hf.integrity(6, 2L, 1L, 2L, 12L, 22L);
    assertTrue(hf.isValid());
    hf.reset();
    hf.integrity(0, 0L, 0L, 0L, 0L, 0L);
  }

  public void testClone() throws CloneNotSupportedException {
    final ProteinMask mask = dummyMask();
    final ProteinMask mask2 = mask.clone();
    assertNotNull(mask2);
    assertTrue(mask2 != mask);
    assertEquals(mask.toString(), mask2.toString());
  }

  public void testInvalidSkel() {
    try {
      ProteinMask.factory(new Skeleton(1, 2, 3, 4, 5));
    } catch (RuntimeException e) {
      assertTrue(e.getMessage().startsWith("Invalid skeleton:"));
    }
  }

  public void testSetValues() {
    try {
      dummyMask().setValues(7, false);
      fail("UnsupportedOperationException expected");
    } catch (UnsupportedOperationException e) {
    }
  }

  public void testFastScore() {
    assertEquals(0, dummyMask().fastScore(7));
  }

  public void testIndelScore() {
    assertEquals(0, dummyMask().indelScore(7));
  }

  public void testReadAllIntLongLong() {
    try {
      dummyMask().readAll(7, 0L, 1L);
      fail("UnsupportedOperationException expected");
    } catch (UnsupportedOperationException e) {
    }
  }

  public void testTemplateAllIntLongLong() {
    try {
      dummyMask().templateAll(7, 0L, 1L);
      fail("UnsupportedOperationException expected");
    } catch (UnsupportedOperationException e) {
    }
  }

  public void testFactory() {
    final Skeleton sk = new Skeleton(4, 4, 0, 0, 1);
    final HashFunctionFactory factory = ProteinMask.factory(sk);
    assertEquals(20, factory.hashBits());
    assertEquals(20, factory.windowBits());
    assertEquals(4, factory.windowSize());
    assertEquals(1, factory.numberWindows());
  }

  private ProteinMask dummyMask() {
    final StringWriter sb = new StringWriter();
    final HashFunctionFactory factory = ProteinMask.factory(new Skeleton(4, 4, 0, 0, 1));
    return (ProteinMask) factory.create(new ReadCallMock(sb), new TemplateCallMock(sb));
  }


  /** protein characters */
  static final char[] CHARS = {'X', '\0', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};

  /**
   * number of protein characters
   * @return number of characters
   */
  public static int charsLength() {
    return CHARS.length;
  }

  /**
   * Get a protein character
   * @param i index of protein character
   * @return the character
   */
  public static char getChar(int i) {
    return CHARS[i];
  }
}
