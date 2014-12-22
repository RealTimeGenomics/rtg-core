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

import java.io.File;
import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;



/**
 */
public class LogFileTest extends TestCase {

  private static final String FS = System.getProperty("file.separator");

  /**
   * Test method for {@link com.rtg.util.io.LogFile}.
   * @throws IOException
   */
  public final void test() throws IOException {
    final File file = File.createTempFile("logFileTest", "log");
    final LogStream ls = new LogFile(file);
    assertTrue(file.exists());
    final String ts = ls.toString();
    assertTrue(ts.startsWith("LogFile "));
    assertTrue(ts.contains(FS + "logFileTest"));
    assertTrue(ts.endsWith("log"));
    ls.stream().println("l1");
    ls.stream().println("l2");
    ls.stream().flush();
    assertEquals("l1" + StringUtils.LS + "l2" + StringUtils.LS, FileUtils.fileToString(file));
    assertEquals(file, ls.file());
    assertEquals(new File(file.getPath()), file);
    ls.removeLog();
    assertTrue(!file.exists());
  }

  /**
   * Test method for {@link com.rtg.util.io.LogFile}.
   * @throws IOException
   */
  public final void testNull() throws IOException {
    final File impossibleLog = FileHelper.createTempDirectory();
    try {
      final LogStream ls = new LogFile(impossibleLog);
      assertNull(ls.stream());
      ls.removeLog();
    } finally {
      FileHelper.deleteAll(impossibleLog);
    }
  }
}

