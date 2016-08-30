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

import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 *
 */
public class SortingSamReaderTest extends TestCase {

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File input = FileHelper.resourceToFile("com/rtg/sam/resources/simplestdout-newheader", new File(dir, "input.sam"));
      try (SortingSamReader reader = SortingSamReader.reader(input, SAMFileHeader.SortOrder.queryname, dir)) {
        final String[] expectedRName = {"0", "0", "1", "1", "2", "2"};
        int i = 0;
        for (SAMRecord rec : reader) {
          assertEquals(expectedRName[i++], rec.getReadName());
        }
      }
    }
  }

}