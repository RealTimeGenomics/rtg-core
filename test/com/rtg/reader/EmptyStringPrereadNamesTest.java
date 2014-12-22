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
package com.rtg.reader;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import junit.framework.TestCase;

/**
 * Test class
 */
public class EmptyStringPrereadNamesTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final EmptyStringPrereadNames thing = new EmptyStringPrereadNames(5);
    assertEquals(5, thing.length());
    assertEquals("", thing.name(3));
    final StringBuilder sb = new StringBuilder();
    thing.writeName(sb, 2);
    assertEquals("", sb.toString());
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    thing.writeName(baos, 0);
    assertEquals("", baos.toString());
    final SimplePrereadNames prni = new SimplePrereadNames();
    prni.setName(0, "");
    prni.setName(1, "");
    prni.setName(2, "");
    prni.setName(3, "");
    prni.setName(4, "");
    assertEquals(prni.calcChecksum(), thing.calcChecksum());
  }
}
