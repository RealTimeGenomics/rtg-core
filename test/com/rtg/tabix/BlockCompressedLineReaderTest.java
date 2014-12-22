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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import com.rtg.util.gzip.GzipUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class BlockCompressedLineReaderTest extends TestCase {

  public void test() throws IOException {
    final File dir = FileUtils.createTempDir("bclr", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz", new File(dir, "readerWindow1.sam.gz"));
      final BlockCompressedLineReader bclr = new BlockCompressedLineReader(new BlockCompressedInputStream(sam));
      try {
        final long firstSeekPos = (44947L << 16) | 22870;
        bclr.seek(firstSeekPos);
        assertEquals(firstSeekPos, bclr.getFilePointer());
        final String line = bclr.readLine();
        assertTrue(line.startsWith("857\t147\tsimulatedSequence2\t32834"));
        assertEquals(firstSeekPos, bclr.getLineFilePointer());
        assertEquals(firstSeekPos + line.length() + 1, bclr.getFilePointer());
        final String line2 = bclr.readLine();
        assertTrue(line2.startsWith("251\t99\tsimulatedSequence2\t33229"));
        assertEquals((int) '9', bclr.peek());
        final String line3 = bclr.readLine();
        assertTrue(line3.startsWith("91\t163\tsimulatedSequence2\t33238"));
        assertEquals(3, bclr.getLineNumber());
      } finally {
        bclr.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testLinearRead() throws IOException {
    final File dir = FileUtils.createTempDir("bclr", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz", new File(dir, "readerWindow1.sam.gz"));
      final BlockCompressedLineReader bclr = new BlockCompressedLineReader(new BlockCompressedInputStream(sam));
      try {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(GzipUtils.createGzipInputStream(new FileInputStream(sam))))) {
          String lineA;
          String lineB;
          while (true) {
            lineA = br.readLine();
            lineB = bclr.readLine();
            if (lineA == null || lineB == null) {
              break;
            }
            assertEquals(lineA, lineB);
          }
          assertNull(lineA);
          assertNull(lineB);
        }
      } finally {
        bclr.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
