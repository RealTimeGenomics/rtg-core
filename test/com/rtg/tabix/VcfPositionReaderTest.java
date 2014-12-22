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
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class VcfPositionReaderTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("test", "vcf");
    try {
      final File vcfFile = BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/vcf.txt", new File(dir, "test.vcf.gz"));
      try (VcfPositionReader reader = new VcfPositionReader(new BlockCompressedLineReader(new BlockCompressedInputStream(vcfFile)), 0)) {
        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(14369, reader.getStartPosition());
        assertEquals(1, reader.getLengthOnReference());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(17329, reader.getStartPosition());
        assertEquals(1, reader.getLengthOnReference());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1110695, reader.getStartPosition());
        assertEquals(1, reader.getLengthOnReference());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1230236, reader.getStartPosition());
        assertEquals(1, reader.getLengthOnReference());

        assertTrue(reader.hasNext());
        reader.next();
        assertEquals("20", reader.getReferenceName());
        assertEquals(1234566, reader.getStartPosition());
        assertEquals(4, reader.getLengthOnReference());

        assertFalse(reader.hasNext());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
