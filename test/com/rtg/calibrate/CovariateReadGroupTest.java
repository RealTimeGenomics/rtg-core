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

package com.rtg.calibrate;

import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateReadGroupTest extends TestCase {

  public void test() {
    final CovariateReadGroup rg = new CovariateReadGroup();
    assertEquals("readgroup", rg.name());
    assertEquals(CovariateEnum.READGROUP, rg.getType());
    assertEquals(0, rg.size());
    assertEquals(0, rg.parse("foo"));
    assertEquals(1, rg.newSize());
    assertEquals(0, rg.size());
    assertEquals(true, rg.sizeChanged());
    rg.resized();
    assertEquals(false, rg.sizeChanged());
    assertEquals(1, rg.size());
    assertEquals(1, rg.parse("bar"));
    assertEquals(2, rg.newSize());
    assertEquals(0, rg.parse("foo"));
    assertEquals(2, rg.newSize());
    assertEquals("foo", rg.valueString(0));
    assertEquals("bar", rg.valueString(1));
    assertEquals(-1, rg.valueOf("dasfdfas"));
  }


  public void testValue() {
    final CovariateReadGroup rg = new CovariateReadGroup();
    final SAMRecord sam = new SAMRecord(null);
    assertEquals(0, rg.value(sam, null));
    assertEquals(1, rg.newSize());
    assertEquals("unspecified", rg.valueString(0));

    sam.setAttribute("RG", "foo");
    assertEquals(1, rg.value(sam, null));
    assertEquals(2, rg.newSize());
    assertEquals("foo", rg.valueString(1));
  }
}
