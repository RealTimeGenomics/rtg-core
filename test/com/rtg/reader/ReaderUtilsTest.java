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
import java.util.Map;
import java.util.UUID;

import com.rtg.mode.SequenceType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * JUnit tests for the ReaderUtils class.
 *
 */
public class ReaderUtilsTest extends TestCase {

  public void testIsPairedEndDirectory() throws IOException {
    final File dir = FileUtils.createTempDir("rut", "dontcare");
    try {
      assertFalse(ReaderUtils.isPairedEndDirectory(null));
      assertFalse(ReaderUtils.isPairedEndDirectory(dir));
      final File left = new File(dir.getAbsolutePath() + File.separator + "left");
      final File right = new File(dir.getAbsolutePath() + File.separator + "right");
      assertTrue(left.mkdir());
      assertFalse(ReaderUtils.isPairedEndDirectory(dir));
      assertTrue(right.mkdir());
      assertTrue(ReaderUtils.isPairedEndDirectory(dir));
      assertEquals(left, ReaderUtils.getLeftEnd(dir));
      assertEquals(right, ReaderUtils.getRightEnd(dir));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testValidateNotEmpty() {
    Diagnostic.setLogStream();
    try {
      ReaderUtils.validateNotEmpty(null);
    } catch (NoTalkbackSlimException e) {
      assertEquals("The SDF \"<Unknown>\" was empty.", e.getMessage());
    }
    try {
      ReaderUtils.validateNotEmpty(new MockSequencesReader(SequenceType.DNA, 0));
    } catch (NoTalkbackSlimException e) {
      assertEquals("The SDF \"<Unknown>\" was empty.", e.getMessage());
    }
    try {
      ReaderUtils.validateNotEmpty(new MockSequencesReader(SequenceType.DNA, 0) {
        @Override
        public File path() {
          return new File("boo");
        }
      });
    } catch (NoTalkbackSlimException e) {
      assertEquals("The SDF \"" + new File("boo").getAbsolutePath() + "\" was empty.", e.getMessage());
    }
    ReaderUtils.validateNotEmpty(new MockSequencesReader(SequenceType.DNA, 1));
  }

  public void testgetGuid() throws IOException {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("readerutils", "getGuid");
    try {
      final File left = new File(dir, "left");
      ReaderTestUtils.getReaderDNA(">A\nACGT", left, new SdfId(5L)).close();
      final File right = new File(dir, "right");
      ReaderTestUtils.getReaderDNA(">A\nACGT", right, new SdfId(5L)).close();

      assertEquals(new SdfId(new UUID(0, 5L)), ReaderUtils.getSdfId(left));
      assertEquals(new SdfId(new UUID(0, 5L)), ReaderUtils.getSdfId(dir));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testSequencNameMap() throws IOException {
    final SequencesReader reader = new MockSequencesReader(SequenceType.DNA, 10, 10);
    checkNameMap(reader);
  }

  private void checkNameMap(SequencesReader reader) throws IOException {
    final Map<String, Long> names = ReaderUtils.getSequenceNameMap(reader);
    assertEquals(10, names.size());
    for (int i = 0; i < names.size(); i++) {
      final String key = "seq" + i;
      assertTrue(names.containsKey(key));
      assertEquals(i, (long) names.get(key));
    }
  }

}
