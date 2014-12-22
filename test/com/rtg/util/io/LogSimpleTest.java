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
package com.rtg.util.io;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;



/**
 */
public class LogSimpleTest extends TestCase {

  /**
   * Test method for {@link com.rtg.util.io.LogSimple}.
   */
  public final void test() {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(ba);
    final LogStream ls = new LogSimple(ps);
    assertEquals("LogSimple", ls.toString());
    ls.stream().println("l1");
    ls.stream().println("l2");
    ls.stream().flush();
    assertEquals("l1" + StringUtils.LS + "l2" + StringUtils.LS, ba.toString());
    ls.removeLog(); //does nothing
    assertNull(ls.file());
  }
}

