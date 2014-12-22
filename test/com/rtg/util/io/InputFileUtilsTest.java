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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

/**
 * JUnit tests for the InputFileUtils class.
 */
public class InputFileUtilsTest extends TestCase {

  public void testRemoveRedundantPaths() throws IOException {
    final List<File> files = new ArrayList<>();
    for (int i = 0; i < 5; i++) {
      final File f = new File("file" + i + ".txt");
      for (int j = 0; j <= i; j++) {
        files.add(f);
      }
    }
    final List<File> cleaned = InputFileUtils.removeRedundantPaths(files);
    assertEquals(15, files.size());
    assertEquals(5, cleaned.size());
    final Iterator<File> it = cleaned.iterator();
    for (int i = 0; i < 5; i++) {
      final File f = it.next();
      assertEquals("file" + i + ".txt", f.getPath());
    }
  }

  public void testCheckIdenticalPaths() throws IOException {
    final File f1 = new File("f1");
    final File f2 = new File("f2");
    assertFalse(InputFileUtils.checkIdenticalPaths(f1, f2));
    assertTrue(InputFileUtils.checkIdenticalPaths(f1, f1));
  }
}
