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
package com.rtg.tabix;

import java.io.File;
import java.io.IOException;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class AlleleCountsPositionReaderTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("test", "vcf");
    try {
      final File acFile = FileHelper.resourceToFile("com/rtg/variant/resources/alleleCounts.ac.gz", new File(dir, "blah.ac.gz"));
      try (AlleleCountsPositionReader reader = new AlleleCountsPositionReader(new BlockCompressedLineReader(new BlockCompressedInputStream(acFile)), 0)) {
        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1, reader.getLengthOnReference());
        assertEquals(14369, reader.getStartPosition());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1, reader.getLengthOnReference());
        assertEquals(17329, reader.getStartPosition());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1, reader.getLengthOnReference());
        assertEquals(1110695, reader.getStartPosition());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1, reader.getLengthOnReference());
        assertEquals(1230236, reader.getStartPosition());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(4, reader.getLengthOnReference());
        assertEquals(1234566, reader.getStartPosition());

        assertFalse(reader.hasNext());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
