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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import net.sf.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class TabixIndexMergeTest extends TestCase {

  private static final String SAM_RESOURCE = "com/rtg/tabix/resources";
  private static final String SAM_FILES = "tabixmerge%d.sam.gz";

  private NanoRegression mNano;
  @Override
  public void setUp() {
    mNano = new NanoRegression(TabixIndexMergeTest.class);
  }
  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }
  public void testSam() throws Exception {
    final File dir = FileUtils.createTempDir("indexmerge", "test");
    try {
      final ArrayList<File> files = new ArrayList<>();
      final ArrayList<Long> dataFileSizes = new ArrayList<>();
      for (int i = 1; i <= 4; i++) {
        final String samFileName = String.format(SAM_FILES, i);
        final File samFile = new File(dir, samFileName);
        final File tbiFile = new File(dir, samFileName + ".tbi");
        FileHelper.resourceToFile(String.format("%s/%s", SAM_RESOURCE, samFileName), samFile);
        FileHelper.resourceToFile(String.format("%s/%s.tbi", SAM_RESOURCE, samFileName), tbiFile);
        files.add(tbiFile);
        dataFileSizes.add(samFile.length());
      }
      final File mergedIndex = new File(dir, "merged.sam.gz.tbi");
      TabixIndexMerge.mergeTabixFiles(mergedIndex, files, dataFileSizes);
      try (InputStream fis = new BlockCompressedInputStream(new FileInputStream(mergedIndex))) {
        final String indexDebug = IndexTestUtils.tbiIndexToUniqueString(fis);
        mNano.check("merged.sam.gz.tbi.debug", indexDebug);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

}
