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

package com.rtg.variant.bayes.multisample.family;

import junit.framework.TestCase;

/**
 */
public class UniqueIdTest extends TestCase {

  public void test() {
    final UniqueId ui = new UniqueId(3);
    ui.globalIntegrity();
    assertEquals(0, ui.numberIdsSoFar());
    assertEquals(-1, ui.id(0));
    assertEquals(-1, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(0, ui.addId(1));
    assertEquals(0, ui.addId(1));
    ui.globalIntegrity();
    assertEquals(1, ui.numberIdsSoFar());
    assertEquals(-1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(1, ui.addId(0));
    assertEquals(1, ui.addId(0));
    ui.globalIntegrity();
    assertEquals(2, ui.numberIdsSoFar());
    assertEquals(1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(2, ui.addId(2));
    assertEquals(2, ui.addId(2));
    ui.globalIntegrity();
    assertEquals(3, ui.numberIdsSoFar());
    assertEquals(1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(2, ui.id(2));

  }
}
