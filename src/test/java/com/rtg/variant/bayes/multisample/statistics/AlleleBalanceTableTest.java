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
