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
import java.io.RandomAccessFile;

import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.SimpleArchive;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 *
 *
 */
public class SequenceStreamManagerTest extends DefaultSequencesReaderTest {
  public SequenceStreamManagerTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(SequenceStreamManagerTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testTruncateFile() throws IOException {
    final File temp = FileUtils.createTempDir("reader", "ssm");
    try {
      try (InputStream is = Resources.getResourceAsStream("com/rtg/reader/resources/sdfver8.arch")) {
        SimpleArchive.unpackArchive(is, mDir);
        final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(mDir);
        try {
          byte[] buff = new byte[(int) sr.maxLength()];
          while (sr.nextSequence()) {
            sr.readCurrent(buff);
          }
        } finally {
          sr.close();
        }
        final RandomAccessFile raf = new RandomAccessFile(new File(mDir, SdfFileUtils.SEQUENCE_DATA_FILENAME + "0"), "rw");
        try {
          raf.setLength(raf.length() - 50);
        } finally {
          raf.close();
        }
      }
      try (SequencesReader sr2 = SequencesReaderFactory.createDefaultSequencesReader(mDir)) {
        byte[] buff2 = new byte[(int) sr2.maxLength()];
        try {
          while (sr2.nextSequence()) {
            sr2.readCurrent(buff2);
          }
          fail("should have dided horribly");
        } catch (final CorruptSdfException exp) {
          //expected
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }
}

