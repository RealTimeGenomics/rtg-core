/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
