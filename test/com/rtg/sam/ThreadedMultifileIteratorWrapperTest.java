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
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;

import htsjdk.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;

/**
 */
public class ThreadedMultifileIteratorWrapperTest extends TestCase {

  public void test() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory top = new TestDirectory()) {
      final File first = new File(top, "first");
      FileUtils.stringToFile(ThreadedMultifileIteratorTest.SAM_HEAD1 + "@RG\tID:1\tSM:sample1\tPL:ILLUMINA" + "\n" + MultifileIteratorTest.SAM_REC_RG1, first);
      final File second = new File(top, "second");
      FileUtils.stringToFile(ThreadedMultifileIteratorTest.SAM_HEAD1 + "@RG\tID:2\tSM:sample1\tPL:ILLUMINA" + "\n" + MultifileIteratorTest.SAM_REC_RG2, second);
      final ArrayList<File> fileList = new ArrayList<>();
      fileList.add(first);
      fileList.add(second);
      final SamFilterParams params = SamFilterParams.builder().create();
      final SingletonPopulatorFactory<VariantAlignmentRecord> pf = new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator());
      try (ThreadedMultifileIteratorWrapper<VariantAlignmentRecord> it = new ThreadedMultifileIteratorWrapper<>(fileList, 1, pf, params, SamUtils.getUberHeader(fileList))) {
        it.setSequenceId(0);
        final List<SAMReadGroupRecord> readGroups = it.header().getReadGroups();
        assertEquals(2, readGroups.size());
        boolean foundFirst = false;
        boolean foundSecond = false;
        for (final SAMReadGroupRecord r : readGroups) {
          if (r.getReadGroupId().equals("1")) {
            foundFirst = true;
          }
          if (r.getReadGroupId().equals("2")) {
            foundSecond = true;
          }
        }
        assertTrue(foundFirst && foundSecond);
        int count = 0;
        while (it.hasNext()) {
          it.next();
          count++;
        }
        assertEquals(2, count);
      }
    }
  }

}
