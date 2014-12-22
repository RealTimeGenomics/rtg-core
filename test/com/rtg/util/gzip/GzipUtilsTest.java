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
package com.rtg.util.gzip;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.zip.GZIPOutputStream;

import com.rtg.util.io.IOUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class GzipUtilsTest extends TestCase {

  public void testTwoStreams() throws IOException {
    final ByteArrayOutputStream firstBaos = new ByteArrayOutputStream();
    final ByteArrayOutputStream secondBaos = new ByteArrayOutputStream();
    final ByteArrayOutputStream joinedBaos = new ByteArrayOutputStream();

    GZIPOutputStream out = new GZIPOutputStream(firstBaos);
    for (byte i = 0; i < 10; i++) {
      out.write(i);
    }
    out.flush();
    out.close();
    out = new GZIPOutputStream(secondBaos);
    for (byte i = 10; i < 20; i++) {
      out.write(i);
    }
    out.flush();
    out.close();
    joinedBaos.write(firstBaos.toByteArray());
    joinedBaos.write(secondBaos.toByteArray());
    joinedBaos.flush();
    joinedBaos.close();


    final InputStream in = GzipUtils.createGzipInputStream(new ByteArrayInputStream(joinedBaos.toByteArray()));
    try {
      for (byte i = 0; i < 20; i++) {
        assertEquals(i, in.read());
      }
      assertEquals(-1, in.read());

    } finally {
      in.close();
    }

    final InputStream in2 = GzipUtils.createGzipInputStream(new ByteArrayInputStream(joinedBaos.toByteArray()));
    try {
      final byte[] b = new byte[3];
      int count = 0;
      int num;
      int i = 0;
      while ((num = in2.read(b, 0, 3)) > 0) {

        count += num;
        for (int j = 0; j < num; j++, i++) {
          assertEquals(i, b[j]);
        }
      }
      assertEquals(20, count);
      assertEquals(20, i);
      final int exp = -1;
      assertEquals(exp,  in2.read(b, 0, 3));
    } finally {
      in.close();
    }
  }

  public void testSingleStream() throws IOException {
    final ByteArrayOutputStream firstBaos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(new GZIPOutputStream(firstBaos));

    ps.print("12345");
    ps.flush();
    ps.close();

    final InputStream in = GzipUtils.createGzipInputStream(new ByteArrayInputStream(firstBaos.toByteArray()));
    try {
      final BufferedReader r = new BufferedReader(new InputStreamReader(in));
      try {
        assertEquals("12345", r.readLine());
        assertNull(r.readLine());
      } finally {
        r.close();
      }
    } finally {
      in.close();
    }
  }

  public void testEmptyStream() throws IOException {
    final ByteArrayOutputStream firstBaos = new ByteArrayOutputStream();
    final GZIPOutputStream out = new GZIPOutputStream(firstBaos);
    out.flush();
    out.close();

    final InputStream in = GzipUtils.createGzipInputStream(new ByteArrayInputStream(firstBaos.toByteArray()));
    try {
      assertTrue(in.read() <= 0);

    } finally {
      in.close();
    }
  }

  public void testShouldOverrideGzip1() throws IOException {
    for (int i = 1; i < GzipUtils.GZIP_FILE.length; i++) {
      final ByteArrayInputStream javaDevelopersAreStupid = new GzipUtils.FakeBlockingInputStream(GzipUtils.GZIP_FILE);
      final byte[] input = IOUtils.readData(GzipUtils.createGzipInputStream(new BufferedInputStream(javaDevelopersAreStupid, i)));
      assertEquals("i: " + i, new String(GzipUtils.EXPECTED_UNZIPPED), new String(input));
    }
  }

  public void testShouldOverrideGzip2() throws IOException {
    final File tmp = FileHelper.createTempDirectory();
    final boolean previousOverride = GzipUtils.getOverrideGzip();
    try {
      final File out = new File(tmp, "out");
      final GZIPOutputStream stream = new GZIPOutputStream(new FileOutputStream(out));
      try {
        final String result = "ch1";
        stream.write(result.getBytes());
        stream.flush();
      } finally {
        stream.close();
      }
      ByteArrayInputStream bais = new ByteArrayInputStream(IOUtils.readData(out));
      GzipUtils.setOverrideGzip(true);
      assertTrue(GzipUtils.createGzipInputStream(bais) instanceof WorkingGzipInputStream);
      bais = new ByteArrayInputStream(IOUtils.readData(out));
      assertTrue(GzipUtils.createGzipInputStream(bais, 1) instanceof WorkingGzipInputStream);

      bais = new ByteArrayInputStream(IOUtils.readData(out));
      GzipUtils.setOverrideGzip(false);
      assertFalse(GzipUtils.createGzipInputStream(bais) instanceof WorkingGzipInputStream);
      bais = new ByteArrayInputStream(IOUtils.readData(out));
      assertFalse(GzipUtils.createGzipInputStream(bais, 1) instanceof WorkingGzipInputStream);
    } finally {
      GzipUtils.setOverrideGzip(previousOverride);
      FileHelper.deleteAll(tmp);
    }
  }

  public void testShouldOverrideGzip3() throws Exception {
    final File tmp = FileHelper.createTempDirectory();
    try {
      final File out = new File(tmp, "out");
      final GZIPOutputStream stream = new GZIPOutputStream(new FileOutputStream(out));
      try {
        final String result = "ch1";
        stream.write(result.getBytes());
        stream.flush();
      } finally {
        stream.close();
      }
      final ByteArrayInputStream bais = new ByteArrayInputStream(IOUtils.readData(out));
      try {
        if (GzipUtils.createGzipInputStream(bais) instanceof WorkingGzipInputStream) {
          assertTrue(GzipUtils.shouldOverrideGzip());
        } else {
          assertFalse(GzipUtils.shouldOverrideGzip());
        }
      } finally {
        bais.close();
      }
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }
}
