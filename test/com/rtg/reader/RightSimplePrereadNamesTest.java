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

import java.io.IOException;

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class RightSimplePrereadNamesTest extends TestCase {
  public void testSomeMethod() throws IOException {
    final SimplePrereadNames sprn = new SimplePrereadNames();
    final RightSimplePrereadNames rprn = new RightSimplePrereadNames(sprn);
    sprn.setName(0, "first");
    sprn.setName(1, "second/1");
    sprn.setName(2, "third");
    sprn.setName(3, "fourth/1");
    rprn.setName(0, "f");
    rprn.setName(1, "sec");
    rprn.setName(2, "thi");
    rprn.setName(3, "fou");
    assertEquals(4L, sprn.length());
    assertEquals(4L, rprn.length());
    assertEquals("first", sprn.name(0));
    assertEquals("second/1", sprn.name(1));
    assertEquals("third", sprn.name(2));
    assertEquals("fourth/1", sprn.name(3));
    assertEquals("first", rprn.name(0));
    assertEquals("second/1", rprn.name(1));
    assertEquals("third", rprn.name(2));
    assertEquals("fourth/1", rprn.name(3));
    assertEquals(66, sprn.bytes());
    assertEquals(0, rprn.bytes());
    StringBuilder sb = new StringBuilder();
    sprn.writeName(sb, 2);
    assertEquals("third", sb.toString());
    MemoryPrintStream mps = new MemoryPrintStream();
    sprn.writeName(mps.outputStream(), 1);
    assertEquals("second/1", mps.toString());

    sb = new StringBuilder();
    rprn.writeName(sb, 2);
    assertEquals("third", sb.toString());
    mps = new MemoryPrintStream();
    rprn.writeName(mps.outputStream(), 1);
    assertEquals("second/1", mps.toString());
}

}
