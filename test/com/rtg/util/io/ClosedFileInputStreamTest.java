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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ClosedFileInputStreamTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("closedStream", "test");
    try {
      final File data = new File(dir, "data");
      final Random r = new Random();
      byte[] hundy = new byte[10 * 1024 * 1024];
      r.nextBytes(hundy);
      FileHelper.streamToFile(new ByteArrayInputStream(hundy), data);
      assertTrue(Arrays.equals(hundy, IOUtils.readData(new ClosedFileInputStream(data))));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void readAndCheckBlock(ClosedFileInputStream cfi, byte[] src, int position) throws IOException {
    final int toread = 100;
    final byte[] expected = new byte[toread];
    System.arraycopy(src, position, expected, 0, toread);

    final byte[] readdata = new byte[toread];
    cfi.seek(position);

    int read = 0;
    int batchread;
    while (read < toread && (batchread = cfi.read(readdata, read, toread - read)) != -1) {
      read += batchread;
    }

    assertTrue(Arrays.equals(expected, readdata));
  }

  public void testSeek() throws IOException {
    final int bufSize = 16 * 1024;
    final int fileSize = 5 * bufSize;
    final File dir = FileUtils.createTempDir("closedStream", "test");
    try {
      final File data = new File(dir, "data");
      final Random r = new Random();
      final byte[] fileContents = new byte[fileSize];
      r.nextBytes(fileContents);
      FileHelper.streamToFile(new ByteArrayInputStream(fileContents), data);

      try (ClosedFileInputStream cfi = new ClosedFileInputStream(data, bufSize)) {
        // First will seek and read from disk
        readAndCheckBlock(cfi, fileContents, 100);
        assertEquals(1, cfi.getDiskSeeksDone());

        // This will seek and read from disk since it's > bufSize away
        readAndCheckBlock(cfi, fileContents, bufSize + 100);
        assertEquals(2, cfi.getDiskSeeksDone());

        // This seek but not read from disk
        readAndCheckBlock(cfi, fileContents, bufSize + 200);
        assertEquals(2, cfi.getDiskSeeksDone());

        // This seek but not read from disk (back to the start read previously)
        readAndCheckBlock(cfi, fileContents, bufSize + 100);
        assertEquals(2, cfi.getDiskSeeksDone());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testSomeMethod2() throws IOException {
    final File dir = FileUtils.createTempDir("closedStream", "test");
    try {
      final File data = new File(dir, "data");
      final String exp = "a line of text is easier to read than line of bytes";
      FileUtils.stringToFile(exp, data);
      assertEquals(exp, IOUtils.readAll(new ClosedFileInputStream(data, 10)));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
