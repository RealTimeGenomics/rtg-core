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
package com.rtg.launcher;

import java.io.File;

import com.rtg.reader.FormatCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class GlobalFlagsTest extends TestCase {

  private static final CFlags FLAGS = new CFlags();

  static final String TEST_FLAG = "test-flag";
  static final Integer TEST_DEFAULT = 20;

  static {
    GlobalFlags.registerFlag(TEST_FLAG, Integer.class, TEST_DEFAULT);
    GlobalFlags.registerExperimentalFlags(FLAGS);
  }

  @Override
  public void setUp() {
    GlobalFlags.resetAccessedStatus();
  }

  public void test() {
    FLAGS.setFlags();
    assertEquals(TEST_DEFAULT, GlobalFlags.getFlag(TEST_FLAG).getValue());

    GlobalFlags.resetAccessedStatus();
    FLAGS.setFlags("--XX" + TEST_FLAG, "902");
    assertEquals(902, GlobalFlags.getFlag(TEST_FLAG).getValue());
    assertEquals(902, GlobalFlags.getIntegerValue(TEST_FLAG));
    assertTrue(GlobalFlags.isSet(TEST_FLAG));
  }

  public void testAccessCheck() {
    assertTrue(GlobalFlags.initialAccessCheck());
    FLAGS.setFlags("--XX" + TEST_FLAG, "902");
    assertTrue(GlobalFlags.initialAccessCheck());
    GlobalFlags.getIntegerValue(TEST_FLAG);
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      assertFalse(GlobalFlags.initialAccessCheck());
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testUnused() {
    FLAGS.setFlags();
    assertTrue(GlobalFlags.finalAccessCheck());
    FLAGS.setFlags("--XX" + TEST_FLAG, "34");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      assertFalse(GlobalFlags.finalAccessCheck());
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testUnusedRealWorld() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final MemoryPrintStream err = new MemoryPrintStream();

    try (TestDirectory td = new TestDirectory()) {
      final File f = new File(td, "f.fq.gz");
      FileHelper.stringToGzFile(">blah\nacgt\n", f);

      final FormatCli cli = new FormatCli();
      cli.mainInit(new String[]{"-o", td.getPath() + "/meh", f.getPath(), "--XX" + TEST_FLAG, "33"}, mps.outputStream(), err.printStream());

      assertTrue(err.toString(), err.toString().contains("XX" + TEST_FLAG + " is set but never accessed"));

    }
  }
}
