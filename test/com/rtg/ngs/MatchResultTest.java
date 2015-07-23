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
package com.rtg.ngs;


import junit.framework.TestCase;

/**
 */
public class MatchResultTest extends TestCase {

  public final void test() {
    final MatchResult r = new MatchResult(0);
    r.addMatchResult(0, 1, 2, true);
    r.addMatchResult(2, 234, 5, true);
    r.addMatchResult(5, 14, 2, true);
    r.addMatchResult(76, 14, 2, true);
    r.addMatchResult(12, 34, 2, false);
    assertEquals(0, r.getTemplateId(0));
    assertEquals(2, r.getTemplateId(1));
    assertEquals(5, r.getTemplateId(2));
    assertEquals(76, r.getTemplateId(3));
    assertEquals(12, r.getTemplateId(4));
    assertEquals(1, r.getPosition(0));
    assertEquals(34, r.getPosition(4));
    assertEquals(2, r.getEncodedReadId(0));
    assertEquals(2, r.getEncodedReadId(4));
    assertEquals(true, r.isReverse(0));
    assertEquals(false, r.isReverse(4));
    r.sort();
    assertEquals(0, r.getTemplateId(0));
    assertEquals(2, r.getTemplateId(1));
    assertEquals(5, r.getTemplateId(2));
    assertEquals(12, r.getTemplateId(3));
    assertEquals(76, r.getTemplateId(4));
    assertEquals(1, r.getPosition(0));
    assertEquals(14, r.getPosition(4));
    assertEquals(2, r.getEncodedReadId(0));
    assertEquals(2, r.getEncodedReadId(4));
    assertEquals(true, r.isReverse(0));
    assertEquals(true, r.isReverse(4));

  }

}
