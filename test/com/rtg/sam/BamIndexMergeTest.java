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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import com.rtg.tabix.IndexTestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class BamIndexMergeTest extends TestCase {

  private NanoRegression mNano;
  @Override
  public void setUp() {
    mNano = new NanoRegression(BamIndexMergeTest.class);
  }
  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private static final String BAM_RESOURCE = "com/rtg/sam/resources";
  private static final String BAM_FILES = "indexmerge%d.bam";

  public void testBam() throws Exception {
    final File dir = FileUtils.createTempDir("indexmerge", "test");
    try {
      final ArrayList<File> files = new ArrayList<>();
      final ArrayList<Long> dataFileSizes = new ArrayList<>();
      for (int i = 1; i <= 4; i++) {
        final String samFileName = String.format(BAM_FILES, i);
        final File bamFile = new File(dir, samFileName);
        final File baiFile = new File(dir, samFileName + ".bai");
        FileHelper.resourceToFile(String.format("%s/%s", BAM_RESOURCE, samFileName), bamFile);
        FileHelper.resourceToFile(String.format("%s/%s.bai", BAM_RESOURCE, samFileName), baiFile);
        files.add(baiFile);
        dataFileSizes.add(bamFile.length());
      }
      final File mergedIndex = new File(dir, "merged.bam.bai");
      BamIndexMerge.mergeBamIndexFiles(mergedIndex, files, dataFileSizes);
      try (InputStream fis = new FileInputStream(mergedIndex)) {
        final String indexDebug = IndexTestUtils.bamIndexToUniqueString(fis);
        mNano.check("merged.bam.bai.debug", indexDebug);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
