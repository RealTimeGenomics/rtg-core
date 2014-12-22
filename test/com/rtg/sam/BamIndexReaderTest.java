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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;

import com.rtg.tabix.VirtualOffsets;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.SAMFileHeader;

import junit.framework.TestCase;

/**
 * Test class
 */
public class BamIndexReaderTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("bamindexreader", "test");
    try {
      final File bam = new File(dir, "mated.bam");
      FileHelper.resourceToFile("com/rtg/sam/resources/mated.bam", bam);
      final File index = new File(dir, "mated.bam.bai");
      FileHelper.resourceToFile("com/rtg/sam/resources/mated.bam.bai", index);
      final SAMFileHeader header = SamUtils.getSingleHeader(bam);
      BamIndexReader tir = new BamIndexReader(index, header.getSequenceDictionary());
      VirtualOffsets positions = tir.getFilePointers(SamRangeUtils.createExplicitReferenceRange(header, new SamRegionRestriction("simulatedSequence", 0, 5000)));
      assertEquals(151L, positions.start(0));
      assertEquals((446375L << 16) | 49187, positions.end(0));
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
}
