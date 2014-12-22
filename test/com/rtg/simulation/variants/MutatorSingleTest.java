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

package com.rtg.simulation.variants;

import com.rtg.util.PortableRandom;
import com.rtg.util.test.NotRandomRandom;

import junit.framework.TestCase;

/**
 */
public class MutatorSingleTest extends TestCase {

  public void testBad() {
    try {
      new MutatorSingle("X@");
      fail();
    } catch (final RuntimeException e) {
      assertEquals("X@", e.getMessage());
    }
  }
  private void expandCheck(String subject, String expected) {
    String actual = MutatorSingle.expandSpec(subject);
    assertEquals(expected, actual);
  }
  public void testExpand() {
    expandCheck("5I3E4D", "IIIIIEEEDDDD");
    expandCheck("5IEIDEE4D", "IIIIIEIDEEDDDD");
    final String unchanged = "IDIEIXIXIXEIEIDIDDIDI";
    expandCheck(unchanged, unchanged);
    try {
      MutatorSingle.expandSpec("5I8");
      fail();
    } catch (RuntimeException e) {
    }
  }

  public void test1() {
    final MutatorSingle ms = new MutatorSingle("XX");
    ms.integrity();
    assertEquals("XX", ms.toString());
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, left.getFirstHaplotype());
    assertEquals("2:CG:CG", left.toString());
    assertEquals("2:TA:TA", right.toString());
  }

  public void test2() {
    final MutatorSingle ms = new MutatorSingle("YY");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, left.getFirstHaplotype());
    assertEquals("2:CG:CG", left.toString());
    assertEquals("2:GT:GT", right.toString());
  }

  public void test3() {
    final MutatorSingle ms = new MutatorSingle("I");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, left.getFirstHaplotype());
    assertEquals("0:A:A", left.toString());
    assertEquals("0:C:C", right.toString());
  }

  public void test4() {
    final MutatorSingle ms = new MutatorSingle("D");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, left.getFirstHaplotype());
    assertEquals("1::", left.toString());
    assertEquals("1::", right.toString());
  }

  public void test5() {
    final MutatorSingle ms = new MutatorSingle("J");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 0, ra, left.getFirstHaplotype());
    assertEquals("0:A:A", left.toString());
    assertEquals("0:G:G", right.toString());
  }

  public void test6() {
    final MutatorSingle ms = new MutatorSingle("==");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1, 2, 3, 4}, 2, ra, null);
    final MutatorResult right = ms.generateMutation(new byte[] {1, 2, 3, 4}, 2, ra, left.getFirstHaplotype());
    assertEquals("2:GT:GT", left.toString());
    assertEquals("2:GT:GT", right.toString());
  }

  public void test1end() {
    final MutatorSingle ms = new MutatorSingle("==");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1}, 0, ra, null);
    assertNull(left);
  }

  public void test1endE() {
    final MutatorSingle ms = new MutatorSingle("EE");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1}, 0, ra, null);
    assertNull(left);
  }

  public void test2end() {
    final MutatorSingle ms = new MutatorSingle("XX");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1}, 0, ra, null);
    assertNull(left);
  }

  public void test3end() {
    final MutatorSingle ms = new MutatorSingle("YY");
    ms.integrity();
    final PortableRandom ra = new NotRandomRandom();
    final MutatorResult left = ms.generateMutation(new byte[] {1}, 0, ra, null);
    assertNull(left);
  }


  public void testRandom() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.random(random));
    assertEquals(2, MutatorSingle.random(random));
    assertEquals(3, MutatorSingle.random(random));
    assertEquals(4, MutatorSingle.random(random));
    assertEquals(1, MutatorSingle.random(random));
  }

  public void testMinusa() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(2, MutatorSingle.minus(random, (byte) 1));
    assertEquals(3, MutatorSingle.minus(random, (byte) 1));
    assertEquals(4, MutatorSingle.minus(random, (byte) 1));
    assertEquals(2, MutatorSingle.minus(random, (byte) 1));
    assertEquals(3, MutatorSingle.minus(random, (byte) 1));
  }

  public void testMinusb() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.minus(random, (byte) 2));
    assertEquals(3, MutatorSingle.minus(random, (byte) 2));
    assertEquals(4, MutatorSingle.minus(random, (byte) 2));
    assertEquals(1, MutatorSingle.minus(random, (byte) 2));
    assertEquals(3, MutatorSingle.minus(random, (byte) 2));
  }

  public void testMinusc() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.minus(random, (byte) 4));
    assertEquals(2, MutatorSingle.minus(random, (byte) 4));
    assertEquals(3, MutatorSingle.minus(random, (byte) 4));
    assertEquals(1, MutatorSingle.minus(random, (byte) 4));
    assertEquals(2, MutatorSingle.minus(random, (byte) 4));
  }

  public void testMinus2a() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(2, MutatorSingle.minus(random, (byte) 1, (byte) 3));
    assertEquals(4, MutatorSingle.minus(random, (byte) 1, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 1, (byte) 3));
    assertEquals(4, MutatorSingle.minus(random, (byte) 1, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 1, (byte) 3));
  }

  public void testMinus2b() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(3, MutatorSingle.minus(random, (byte) 1, (byte) 2));
    assertEquals(4, MutatorSingle.minus(random, (byte) 1, (byte) 2));
    assertEquals(3, MutatorSingle.minus(random, (byte) 1, (byte) 2));
    assertEquals(4, MutatorSingle.minus(random, (byte) 1, (byte) 2));
    assertEquals(3, MutatorSingle.minus(random, (byte) 1, (byte) 2));
  }

  public void testMinus2c() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.minus(random, (byte) 3, (byte) 4));
    assertEquals(2, MutatorSingle.minus(random, (byte) 3, (byte) 4));
    assertEquals(1, MutatorSingle.minus(random, (byte) 3, (byte) 4));
    assertEquals(2, MutatorSingle.minus(random, (byte) 3, (byte) 4));
    assertEquals(1, MutatorSingle.minus(random, (byte) 3, (byte) 4));
  }

  public void testMinus2d() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.minus(random, (byte) 3, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 3, (byte) 3));
    assertEquals(4, MutatorSingle.minus(random, (byte) 3, (byte) 3));
    assertEquals(1, MutatorSingle.minus(random, (byte) 3, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 3, (byte) 3));
  }

  public void testMinus2e() {
    final PortableRandom random = new NotRandomRandom();
    assertEquals(1, MutatorSingle.minus(random, (byte) 4, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 4, (byte) 3));
    assertEquals(1, MutatorSingle.minus(random, (byte) 4, (byte) 3));
    assertEquals(2, MutatorSingle.minus(random, (byte) 4, (byte) 3));
    assertEquals(1, MutatorSingle.minus(random, (byte) 4, (byte) 3));
  }


}
