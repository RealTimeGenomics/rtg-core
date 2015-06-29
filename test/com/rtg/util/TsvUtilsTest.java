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
