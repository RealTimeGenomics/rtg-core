
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

package com.rtg.util.io.bzip2;

import java.io.IOException;

import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * test class
 */
public class CBZip2InputStreamTest extends TestCase {

  public void testSomeMethod() throws IOException {
    try (CBZip2InputStream s = new CBZip2InputStream(Resources.getResourceAsStream("com/rtg/util/io/bzip2/resources/textfile.bz2"))) {
      final String bzString = FileUtils.streamToString(s);
      final String expString = FileHelper.resourceToString("com/rtg/util/io/bzip2/resources/textfile");
      assertEquals(expString, bzString);
    }
  }

  public void testMultimember() throws IOException {
    try (CBZip2InputStream s = new CBZip2InputStream(Resources.getResourceAsStream("com/rtg/util/io/bzip2/resources/textfilemulti.bz2"))) {
      final String bzString = FileUtils.streamToString(s);
      final String expString = FileHelper.resourceToString("com/rtg/util/io/bzip2/resources/textfile");
      assertEquals(expString, bzString);
    }
  }

  public void testBadFile() throws IOException {
    try (CBZip2InputStream s = new CBZip2InputStream(Resources.getResourceAsStream("com/rtg/util/io/bzip2/resources/textfilebad.bz2"))) {
      try {
        FileUtils.streamToString(s);
        fail();
      } catch (final IOException e) {
        assertTrue(e.getMessage().contains("crc"));
        // expected
      }
    }
  }

  public void testRepetitiveFile() throws IOException {
    try (CBZip2InputStream s = new CBZip2InputStream(Resources.getResourceAsStream("com/rtg/util/io/bzip2/resources/sample3.ref.bz2"))) {
      final String bzString = FileUtils.streamToString(s);
      final String expString = FileHelper.resourceToString("com/rtg/util/io/bzip2/resources/sample3.ref");
      assertEquals(expString, bzString);
    }
  }
}
