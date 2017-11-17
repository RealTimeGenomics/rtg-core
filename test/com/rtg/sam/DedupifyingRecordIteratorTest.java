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
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;

import junit.framework.TestCase;

/**
 * Test class
 */
public class DedupifyingRecordIteratorTest extends TestCase {
  private static final String SAM_RESOURCE = "com/rtg/variant/cnv/resources/testFilter.sam";
  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory("dedup")) {
      final File sam = FileHelper.resourceToFile(SAM_RESOURCE, new File(dir, "testFilter.sam"));
      final List<File> f = new ArrayList<>();
      f.add(sam);
      final SingletonPopulatorFactory<VariantAlignmentRecord> pf = new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0));
      final SamReadingContext context = new SamReadingContext(f, 1, SamFilterParams.builder().excludeUnmapped(true).create(), SamUtils.getUberHeader(f), null);
      try (ThreadedMultifileIterator<VariantAlignmentRecord> tmfi = new ThreadedMultifileIterator<>(context, pf)) {
        final DedupifyingRecordIterator<VariantAlignmentRecord> dd = new DedupifyingRecordIterator<>(tmfi);
        while (dd.hasNext()) {
          dd.next();
          assertEquals(tmfi.getFilteredRecordsCount(), dd.getFilteredRecordsCount());
          assertEquals(tmfi.getInvalidRecordsCount(), dd.getInvalidRecordsCount());
          assertEquals(tmfi.getOutputRecordsCount(), dd.getOutputRecordsCount());
          assertEquals(tmfi.header(), dd.header());
        }
        assertEquals(tmfi.getFilteredRecordsCount(), dd.getFilteredRecordsCount());
        assertEquals(tmfi.getInvalidRecordsCount(), dd.getInvalidRecordsCount());
        assertEquals(tmfi.getOutputRecordsCount(), dd.getOutputRecordsCount());
        assertEquals(tmfi.header(), dd.header());
        dd.close();
        assertEquals(tmfi.getFilteredRecordsCount(), dd.getFilteredRecordsCount());
        assertEquals(tmfi.getInvalidRecordsCount(), dd.getInvalidRecordsCount());
        assertEquals(tmfi.getOutputRecordsCount(), dd.getOutputRecordsCount());
        assertEquals(tmfi.header(), dd.header());
      }
    }
  }
}
