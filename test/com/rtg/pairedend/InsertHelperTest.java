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
package com.rtg.pairedend;

import junit.framework.TestCase;

/**
 */
public class InsertHelperTest extends TestCase {

  /**
   */
  public InsertHelperTest(String name) {
    super(name);
  }

  //  public final void testCalculateInsertSize() {
  //    int insert = InsertHelper.calculateInsertSize(1, false, 5, 10, false, 5);
  //    assertEquals(9, insert);
  //    insert = InsertHelper.calculateInsertSize(1, false, 5, 10, true, 5);
  //    assertEquals(14, insert);
  //    insert = InsertHelper.calculateInsertSize(1, true, 5, 10, false, 5);
  //    assertEquals(4, insert);
  //    insert = InsertHelper.calculateInsertSize(1, true, 5, 10, true, 5);
  //    assertEquals(9, insert);
  //
  //    insert = InsertHelper.calculateInsertSize(10, false, 5, 1, false, 5);
  //    assertEquals(9, insert);
  //    insert = InsertHelper.calculateInsertSize(10, false, 5, 1, true, 5);
  //    assertEquals(4, insert);
  //    insert = InsertHelper.calculateInsertSize(10, true, 5, 1, false, 5);
  //    assertEquals(14, insert);
  //    insert = InsertHelper.calculateInsertSize(10, true, 5, 1, true, 5);
  //    assertEquals(9, insert);
  //  }

  public final void testCalculateTemplateLength() {
    //                                                                      01234567890123456789
    int tlen = InsertHelper.calculateFragmentLength(0, 4, 11, 4);   //  ACGT       ACGT
    assertEquals(15, tlen);
    tlen = InsertHelper.calculateFragmentLength(1, 4, 12, 4);   //       ACGT       ACGT
    assertEquals(15, tlen);

    tlen = InsertHelper.calculateFragmentLength(12, 4, 0, 4);
    assertEquals(16, tlen);
    tlen = InsertHelper.calculateFragmentLength(13, 4, 1, 4);
    assertEquals(16, tlen);
  }

  public void testTlen() {
    checkTlen(150, true, 1, 100, 51, 100);
    checkTlen(100, true, 1, 100, 1, 100);
  }

  private void checkTlen(final int exp, boolean first1, int templateStart1, int readLen1, int templateStart2, int readLen2) {
    checkTlen1(exp, first1, templateStart1, readLen1, templateStart2, readLen2);
    checkTlen1(-exp, !first1, templateStart2, readLen2, templateStart1, readLen1);
  }

  private void checkTlen1(final int exp, boolean first1, int templateStart1, int readLen1, int templateStart2, int readLen2) {
    final int tlen = InsertHelper.tlen(first1, templateStart1, readLen1, templateStart2, readLen2);
    assertEquals(exp, tlen);
  }
}
