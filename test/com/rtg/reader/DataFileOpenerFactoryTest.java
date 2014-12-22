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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import com.rtg.mode.SequenceType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.SeekableStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class DataFileOpenerFactoryTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("temp", "temp");
    try {
      ReaderTestUtils.getDNADir(">abcde\nnacgt", dir);
      final DataFileOpenerFactory fact = new DataFileOpenerFactory(IndexFile.SEQUENCE_ENCODING_COMPRESSED, IndexFile.QUALITY_ENCODING_COMPRESSED, SequenceType.DNA);
      final byte[] b = new byte[5];
      try (InputStream is = fact.getLabelOpener().open(new File(dir, "namedata0"), 6)) {
        final int len = is.read(b);
        assertEquals("abcde".substring(0, len), new String(b, 0, len));
      }
      try (SeekableStream s = fact.getSequenceOpener().openRandomAccess(new File(dir, "seqdata0"), 5)) {
        s.seek(1);
        IOUtils.readFully(s, b, 0, 4);
        assertTrue(Arrays.equals(new byte[]{1, 2, 3, 4, (byte) 'e'}, b));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testOtherMethods() throws Exception {
    final File dir = FileUtils.createTempDir("temp", "temp");
    try {
      final File blah = new File(dir, "blah");
      assertTrue(blah.createNewFile());
      final DataFileOpenerFactory fact = new DataFileOpenerFactory(IndexFile.SEQUENCE_ENCODING_COMPRESSED, IndexFile.QUALITY_ENCODING_COMPRESSED, SequenceType.DNA);
      InputStream is = fact.getSequenceOpener().open(blah, 0);
      try {
        assertTrue(fact.getSequenceOpener().getClass().getName(), is instanceof FileBitwiseInputStream);
      } finally {
        is.close();
      }
      is = fact.getSequenceOpener().openRandomAccess(blah, 0);
      try {
        assertTrue(fact.getSequenceOpener().getClass().getName(), is instanceof FileBitwiseInputStream);
      } finally {
        is.close();
      }
      is = fact.getQualityOpener().open(blah, 0);
      try {
        assertTrue(fact.getQualityOpener().getClass().getName(), is instanceof FileCompressedInputStream);
      } finally {
        is.close();
      }
      is = fact.getQualityOpener().openRandomAccess(blah, 0);
      try {
        assertTrue(fact.getQualityOpener().getClass().getName(), is instanceof FileCompressedInputStream);
      } finally {
        is.close();
      }
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

}
