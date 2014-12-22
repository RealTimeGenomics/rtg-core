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

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.SequenceType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test the corresponding class.
 */
public class TsvSequenceDataSourceTest extends TestCase {

  public static Test suite() {
    return new TestSuite(TsvSequenceDataSourceTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  private File mSample = null;
  private File mSampleBzip2 = null;
  private File mDir = null;
  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileUtils.createTempDir("tsvseqdatasource", "Test");
    mSample = new File(mDir, "sample.tsv.gz");
    FileHelper.stringToGzFile(FileHelper.resourceToString("com/rtg/reader/resources/sample.tsv"), mSample);
    mSampleBzip2 = new File(mDir, "sample.tsv.bz2");
    FileHelper.resourceToFile("com/rtg/reader/resources/sample.tsv.bz2", mSampleBzip2);
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mDir);
    mSample = null;
    mSampleBzip2 = null;
    mDir = null;
  }

  private void checkRead(TsvSequenceDataSource ds, String expected) {
    assertEquals(expected.length(), ds.currentLength());
    final DNAFastaSymbolTable t = new DNAFastaSymbolTable();
    final byte[] b = ds.sequenceData();
    for (int i = 0; i < expected.length(); i++) {
      assertEquals(t.scanResidue(expected.charAt(i)).ordinal(), b[i]);
    }
  }

  private void checkQuality(TsvSequenceDataSource ds, String expected) {
    assertEquals(expected.length(), ds.currentLength());
    final byte[] q = ds.qualityData();
    assertNotNull(q);
    for (int i = 0; i < expected.length(); i++) {
      assertEquals(expected.charAt(i) - '!', q[i]);
    }
  }

  public void test1() throws Exception {
    try (TsvSequenceDataSource ds = new TsvSequenceDataSource(mSample, null)) {
      assertEquals(SequenceType.DNA, ds.type());
      assertTrue(ds.hasQualityData());
      assertTrue(ds.nextSequence());
      assertEquals("1-A", ds.name());
      checkRead(ds, "AAAAAAAAAA");
      checkQuality(ds, "ABCDEFGHIJ");
      assertTrue(ds.nextSequence());
      checkRead(ds, "CCCCCCCCCC");
      checkQuality(ds, "JIHGFEDCBA");
      assertEquals("1-B", ds.name());
      assertTrue(ds.nextSequence());
      checkRead(ds, "GGGGGGGGGG");
      assertEquals("2-A", ds.name());
      assertTrue(ds.nextSequence());
      checkRead(ds, "TTTTTTTTTT");
      assertEquals("2-B", ds.name());
      assertTrue(ds.nextSequence());
      checkQuality(ds, "56789:;<=>");
      checkRead(ds, "ATATATATAT");
      assertTrue(ds.nextSequence());
      checkQuality(ds, "?@ABCDEFGH");
      checkRead(ds, "ATATATATAT");
      assertFalse(ds.nextSequence());
    }
  }

  public void testBzip2() throws Exception {
    try (TsvSequenceDataSource ds = new TsvSequenceDataSource(mSampleBzip2, null)) {
      assertEquals(SequenceType.DNA, ds.type());
      assertTrue(ds.hasQualityData());
      assertTrue(ds.nextSequence());
      assertEquals("1-A", ds.name());
      checkRead(ds, "AAAAAAAAAA");
      checkQuality(ds, "ABCDEFGHIJ");
      assertTrue(ds.nextSequence());
      checkRead(ds, "CCCCCCCCCC");
      checkQuality(ds, "JIHGFEDCBA");
      assertEquals("1-B", ds.name());
      assertTrue(ds.nextSequence());
      checkRead(ds, "GGGGGGGGGG");
      assertEquals("2-A", ds.name());
      assertTrue(ds.nextSequence());
      checkRead(ds, "TTTTTTTTTT");
      assertEquals("2-B", ds.name());
      assertTrue(ds.nextSequence());
      checkQuality(ds, "56789:;<=>");
      checkRead(ds, "ATATATATAT");
      assertTrue(ds.nextSequence());
      checkQuality(ds, "?@ABCDEFGH");
      checkRead(ds, "ATATATATAT");
      assertFalse(ds.nextSequence());
    }
  }

}
