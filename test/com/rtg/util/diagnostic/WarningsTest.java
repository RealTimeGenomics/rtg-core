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

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class WarningsTest extends TestCase {

  public void test() {
    Diagnostic.clearListeners();
    final MemoryPrintStream err = new MemoryPrintStream();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(err.printStream(), out.printStream());
    Diagnostic.addListener(listener);
    Diagnostic.setLogStream();

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());

    final Warnings warnings = new Warnings();
    final Warnings.Warning spy1 = warnings.create(2, "foo1", false);
    final Warnings.Warning spy2 = warnings.create(1, "bar2", true);

    spy1.warn("A");
    spy2.warn("B");
    spy2.warn("C");

    warnings.report();
    mps.close();

    final String t = mps.toString().trim();
    //System.err.println(t);
    assertTrue(t.contains("foo1 A"));
    assertTrue(t.contains("bar2 B"));
    assertTrue(t.contains("bar2 C"));
    assertTrue(t.contains("bar2 occurred 2 times"));

    err.close();
    final String e = err.toString();
    //System.err.println(e);
    assertTrue(e.contains("foo1 A"));
    assertTrue(e.contains("bar2 B"));
    assertFalse(e.contains("bar2 C"));
    assertTrue(e.contains("bar2 occurred 2 times"));

    out.close();
    assertEquals("", out.toString());
  }

}
