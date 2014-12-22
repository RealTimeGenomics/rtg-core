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
package com.rtg.util.diagnostic;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 *
 */
public class NoTalkbackSlimExceptionTest extends TestCase {

  public NoTalkbackSlimExceptionTest(String name) {
    super(name);
  }

  public final void testNoTalkbackSlimExceptionString() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream perr = new PrintStream(bos);
    Diagnostic.setLogStream(perr);
    try {
      final NoTalkbackSlimException e = new NoTalkbackSlimException(new RuntimeException(), ErrorType.INFO_ERROR, "really really long message");
      e.logException();
      perr.flush();
      //System.err.println(bos.toString());
      TestUtils.containsAll(bos.toString(), "really really long message");
    } finally {
      Diagnostic.setLogStream();
      perr.close();
    }

  }

  public final void testNoTalkbackSlimExceptionString2() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream perr = new PrintStream(bos);
    Diagnostic.setLogStream(perr);
    try {
      final NoTalkbackSlimException e = new NoTalkbackSlimException("really really long message");
      e.logException();
      perr.flush();
      //System.err.println(bos.toString());
      TestUtils.containsAll(bos.toString(), "really really long message");

    } finally {
      Diagnostic.setLogStream();
      perr.close();
    }

  }
}
