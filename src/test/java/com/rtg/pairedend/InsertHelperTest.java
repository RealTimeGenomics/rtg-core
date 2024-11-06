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
