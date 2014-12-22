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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.BuildTestUtils;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SimplePrereadNamesTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final SimplePrereadNames sprn = new SimplePrereadNames();
    sprn.setName(0, "first");
    sprn.setName(1, "second");
    sprn.setName(2, "third");
    sprn.setName(3, "fourth");
    assertEquals(4L, sprn.length());
    assertEquals("first", sprn.name(0));
    assertEquals("second", sprn.name(1));
    assertEquals("third", sprn.name(2));
    assertEquals("fourth", sprn.name(3));
    assertEquals(62, sprn.bytes());
    StringBuilder sb = new StringBuilder();
    sprn.writeName(sb, 2);
    assertEquals("third", sb.toString());
    MemoryPrintStream mps = new MemoryPrintStream();
    sprn.writeName(mps.outputStream(), 1);
    assertEquals("second", mps.toString());
  }

  public void testViaFile() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileHelper.createTempDirectory();
    try {
      final File queryDir = BuildTestUtils.prereadDNA(dir, PrereadNamesTest.SEQ_DNA_A2);
      final PrereadNames names = new PrereadNames(queryDir, LongRange.NONE);
      final SimplePrereadNames sprn = new SimplePrereadNames();
      for (long i = 0; i < names.length(); i++) {
        sprn.setName(i, names.name(i));
      }
      assertEquals(names.calcChecksum(), sprn.calcChecksum());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }


}
