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

import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test for corresponding class
 */
public class SimpleArchiveTest extends TestCase {

  private static final String CONTENTS_1 = "a small brown fox";
  private static final String FILENAME_1 = "testFile1";
  private static final String CONTENTS_2 = "lived in a hole";
  private static final String FILENAME_2 = "testFile2";

  public void test() throws IOException {
    final File testDir = FileUtils.createTempDir("simplearchive", "test");
    try {
      final File orig = new File(testDir, "orig");
      final File archive = new File(testDir, "test.dwa");
      final File origFile1 = new File(orig, FILENAME_1);
      final File origFile2 = new File(orig, FILENAME_2);
      assertTrue(orig.mkdir());
      FileUtils.stringToFile(CONTENTS_1, origFile1);
      FileUtils.stringToFile(CONTENTS_2, origFile2);
      SimpleArchive.writeArchive(archive, origFile1, origFile2);
      final File extractDir = new File(testDir, "extract");
      SimpleArchive.unpackArchive(archive, extractDir);
      assertEquals(CONTENTS_1, FileUtils.fileToString(new File(extractDir, FILENAME_1)));
      assertEquals(CONTENTS_2, FileUtils.fileToString(new File(extractDir, FILENAME_2)));
    } finally {
      assertTrue(FileHelper.deleteAll(testDir));
    }
  }


}
