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
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.variant.DefaultMachineErrorChooser;
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
      final File second = new File(top, "second.sam");
      FileUtils.stringToFile(ThreadedMultifileIteratorTest.SAM_HEAD1 + "@RG\tID:2\tSM:sample1\tPL:ILLUMINA" + "\n" + MultifileIteratorTest.SAM_REC_RG2, second);
      final ArrayList<File> fileList = new ArrayList<>();
      fileList.add(first);
      fileList.add(second);
      final SamFilterParams params = SamFilterParams.builder().create();
      final SingletonPopulatorFactory<VariantAlignmentRecord> pf = new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0));
      try (ThreadedMultifileIteratorWrapper<VariantAlignmentRecord> it = new ThreadedMultifileIteratorWrapper<>(new SamReadingContext(fileList, 1, params, SamUtils.getUberHeader(fileList), null), pf)) {
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
          ++count;
        }
        assertEquals(2, count);
      }
    }
  }

}
