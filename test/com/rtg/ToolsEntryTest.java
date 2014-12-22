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
package com.rtg;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class ToolsEntryTest extends TestCase {

  public void test() {
    assertTrue(new ToolsEntry().getSlimModule("HELP") instanceof Command);
  }

  public void test2() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    new ToolsEntry().help(null, mps.outputStream(), mps.printStream());
    final String out = mps.toString();
    TestUtils.containsAll(out, "format", "index", "bgzip");
    assertFalse(out.contains("map"));
    assertFalse(out.contains("family"));
  }
}
