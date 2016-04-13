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

package com.rtg.blacklist;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.AbstractTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

public class HashDistTest extends AbstractTest {

  public void testHistogram() throws IOException {
    final File dir = FileUtils.createTempDir("duster", "test");
    try {
      final File f = new File(dir, "f");
      ReaderTestUtils.getReaderDNA(">t\ntttttttttt", f, null).close();
      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File foo = new File(dir, "foo");
      final int code = new HashDistCli().mainInit(new String[] {f.getPath(), "-o", foo.getPath(), "-s", "1", "-w", "1", "--blacklist-threshold", "1"}, out, err.printStream());
      assertEquals(err.toString(), 0, code);
      err.close();
      assertEquals("10 1" + StringUtils.LS, FileUtils.fileToString(new File(foo, "histogram.txt")));
      assertEquals("T\t10" + StringUtils.LS, FileUtils.fileToString(new File(foo, "blacklist")));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String EXP_HISTO2 = ""
          + "6 1" + StringUtils.LS
          + "7 1" + StringUtils.LS
          + "8 1" + StringUtils.LS
          + "9 1" + StringUtils.LS;

  public void testHistogram2() throws IOException {
    final File dir = FileUtils.createTempDir("duster", "test");
    try {
      final File f = new File(dir, "f");
      //9 as, 8 cs, 6 gs, 7 ts
      ReaderTestUtils.getReaderDNA(">t\nacgtacgtacgtacagtcatcgatcgatca", f, null).close();
      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File foo = new File(dir, "foo");
      final int code = new HashDistCli().mainInit(new String[] {f.getPath(), "-o", foo.getPath(), "-s", "1", "-w", "1", "--blacklist-threshold", "1"}, out, err.printStream());
      assertEquals(err.toString(), 0, code);
      err.close();
      assertEquals(EXP_HISTO2, FileUtils.fileToString(new File(foo, "histogram.txt")));
      TestUtils.containsAll(FileUtils.fileToString(new File(foo, "blacklist")),
        "A\t9" + StringUtils.LS,
        "G\t6" + StringUtils.LS,
        "C\t8" + StringUtils.LS,
        "T\t7" + StringUtils.LS);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
