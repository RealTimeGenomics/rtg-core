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

package com.rtg.variant.bayes.multisample.statistics;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class AlleleBalanceTableTest extends TestCase {

  public void test() {
    final AlleleBalanceTable table = new AlleleBalanceTable();
    table.globalIntegrity();
    assertEquals(0, table.length());

    table.update(0, 0);
    table.update(0, 1);
    table.update(0, 1);
    table.update(1, 1);
    table.update(0, 2);
    table.update(2, 2);
    table.update(2, 2);
    table.update(0, 3);
    table.globalIntegrity();
    assertEquals(7, table.length());

    final String exp = ""
        + "#Coverage R X Count Total" + LS
        + "0 0 0 1 1" + LS
        + "1 0 1 2 3" + LS
        + "1 1 0 1 3" + LS
        + "2 0 2 1 3" + LS
        + "2 2 0 2 3" + LS
        + "3 0 3 1 1" + LS
        ;
    assertEquals(exp, table.toString());
  }

  public void testEmpty() {
    final AlleleBalanceTable table = new AlleleBalanceTable();
    table.globalIntegrity();

    final String exp = ""
        + "#Coverage R X Count Total" + LS
        ;
    assertEquals(exp, table.toString());
  }
}
