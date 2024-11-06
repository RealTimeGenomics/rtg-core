/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
