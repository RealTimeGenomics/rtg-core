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

import java.util.Arrays;

import com.rtg.util.integrity.Exam.ExamException;

import junit.framework.TestCase;

/**
 */
public class MutatorResultTest extends TestCase {

  public void test0() {
    final MutatorResult mr = new MutatorResult(new byte[] {}, new byte[] {}, 0);
    mr.integrity();
    assertEquals("0::", mr.toString());
  }

  public void test1() {
    final MutatorResult mr = new MutatorResult(new byte[] {0, 1, 2}, new byte[] {3, 4}, 2);
    mr.integrity();
    assertEquals("2:NAC:GT", mr.toString());
    assertEquals(2, mr.getConsumed());
    assertTrue(Arrays.equals(new byte[] {0, 1, 2}, mr.getFirstHaplotype()));
    assertTrue(Arrays.equals(new byte[] {3, 4}, mr.getSecondHaplotype()));
  }

  public void testCheckHaplotype() {
    failHaplotype(new byte[]{-1});
    failHaplotype(new byte[]{5});
  }

  private void failHaplotype(byte[] haplotypes) {
    try {
      MutatorResult.checkHaplotype(haplotypes);
      fail();
    } catch (final ExamException e) {
      // expected
    }
  }
}
