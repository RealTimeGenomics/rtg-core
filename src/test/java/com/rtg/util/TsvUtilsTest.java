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
package com.rtg.util;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class TsvUtilsTest extends TestCase {

  private static final String TXT = ""
    + "0" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "155" + TAB + "6M1D1M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACG" + TAB + "```````" + TAB + "AS:i:0" + LS
    + "1" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "155" + TAB + "6M1D1M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACG" + TAB + "```````" + TAB + "AS:i:0" + LS
    + "2" + TAB + "0" + TAB + "g1" + TAB +  "5" + TAB + "155" + TAB + "3M1D4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTC" + TAB + "```````" + TAB + "AS:i:1" + LS
    + "4" + TAB + "0" + TAB + "g1" + TAB +  "6" + TAB + "155" + TAB + "3M1D4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTC" + TAB + "```````" + TAB + "AS:i:1" + LS
    + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "155" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + LS
    ;
  private static final String HEADER = "" + "#a" + TAB + "b" + TAB + "b" + LS;

  public void testCatFiles() throws IOException {
    final File main = FileUtils.createTempDir("tsvtests", "cat");
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final File one = new File(main, "0.tsv");
        final File two = new File(main, "1.tsv");
        FileUtils.stringToFile(HEADER + TXT, one);
        FileUtils.stringToFile(HEADER + TXT, two);
        TsvUtils.tsvCat(bos, one, two);
        bos.flush();
      } finally {
        bos.close();
      }
      assertEquals(HEADER + TXT + TXT, bos.toString());

    } finally {
      assertTrue(FileHelper.deleteAll(main));
    }
  }

  public void testCatFilesSometimesHeader() throws IOException { // Check also works when some files don't have headers
    final File main = FileUtils.createTempDir("tsvtests", "cat");
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final File one = new File(main, "0.tsv");
        final File two = new File(main, "1.tsv");
        final File three = new File(main, "2.tsv");
        FileUtils.stringToFile(HEADER + TXT, one);
        FileUtils.stringToFile(TXT, two);
        FileUtils.stringToFile(TXT, three);
        TsvUtils.tsvCat(bos, one, two, three);
        bos.flush();
      } finally {
        bos.close();
      }
      assertEquals(HEADER + TXT + TXT + TXT, bos.toString());

    } finally {
      assertTrue(FileHelper.deleteAll(main));
    }
  }

  public void testCatFilesGZ() throws IOException {
    final File main = FileUtils.createTempDir("tsvtests", "cat");
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final File one = new File(main, "0.gz");
        final File two = new File(main, "1.gz");
        FileHelper.stringToGzFile(HEADER + TXT, one);
        FileHelper.stringToGzFile(HEADER + TXT, two);
        TsvUtils.tsvCat(bos, one, two);
        bos.flush();
      } finally {
        bos.close();
      }
      assertEquals(HEADER + TXT + TXT, bos.toString());

    } finally {
      assertTrue(FileHelper.deleteAll(main));
    }
  }

}
