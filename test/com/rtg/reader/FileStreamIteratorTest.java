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

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class FileStreamIteratorTest extends TestCase {

  public void testIterator() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("filestreamiterator", "test");
    final int numFiles = 10;
    try {
      final ArrayList<File> files = new ArrayList<>();
      for (int i = 0; i < numFiles; i++) {
        files.add(new File(tempDir, String.valueOf(i)));
        if (files.get(i).createNewFile()) {
          FileUtils.appendToFile(files.get(i), String.valueOf(i));
        }
      }
      final FileStreamIterator iterator = new FileStreamIterator(files, PrereadArm.UNKNOWN);
      try {
        iterator.remove();
        fail();
      } catch (UnsupportedOperationException e) {
        //Expected
      }
      assertNull(iterator.currentFile());
      for (int i = 0; i < numFiles; i++) {
        final InputStream current = iterator.next();
        assertNotNull(current);
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(current))) {
          assertEquals(String.valueOf(i), reader.readLine());
          if (i < numFiles - 1) {
            assertTrue(iterator.hasNext());
            assertTrue(iterator.hasNext());
          } else {
            assertFalse(iterator.hasNext());
          }
          assertEquals(files.get(i), iterator.currentFile());
        }
      }
      InputStream is;
      do {
        is = iterator.next();
      } while (is != null);
      assertNull(iterator.next());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testErrors() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(bos);
    Diagnostic.setLogStream(ps);
    final File tempDir = FileUtils.createTempDir("filestreamiterator", "test");
    try {
      final ArrayList<File> files = new ArrayList<>();
      files.add(new File(tempDir, "nothere"));
      files.add(new File(tempDir, "fake.gz"));
      if (files.get(1).createNewFile()) {
        FileUtils.appendToFile(files.get(1), "blah");
      }
      files.add(new File(tempDir, "real"));
      if (files.get(2).createNewFile()) {
        FileUtils.appendToFile(files.get(2), "blah");
      }
      final FileStreamIterator iterator = new FileStreamIterator(files, null);
      try {
        iterator.hasNext();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage(), e.getMessage().contains("The file: \"" + files.get(0).getPath() + "\" either could not be found or could not be opened."));
      }
      try {
        iterator.next();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage(), e.getMessage().contains("The file: \"" + files.get(1).getPath() + "\" had a problem while reading."));
      }
      iterator.next().close();
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
      Diagnostic.setLogStream();
    }
  }

}
