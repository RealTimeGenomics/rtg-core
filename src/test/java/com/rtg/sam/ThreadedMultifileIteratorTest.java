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

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.Populator;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;

/**
 */
public class ThreadedMultifileIteratorTest extends MultifileIteratorTest {
  private File mDir = null;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
  }
  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  @Override
  RecordIterator<SAMRecord> getIterator(Collection<File> files) throws IOException {
    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    return new ThreadedMultifileIterator<>(new SamReadingContext(files, 2, SamFilterParams.builder().create(), SamUtils.getUberHeader(files), null), pf);
  }

  public static RecordIterator<SAMRecord> getIterator(File dir, String resource) throws IOException {
    return getIterator(dir, SamFilterParams.builder().create(), resource, false);
  }

  public static RecordIterator<SAMRecord> getIterator(File dir, SamFilterParams filterParams, String resource, boolean indexed) throws IOException {
    final File samFile = new File(dir, resource.substring(resource.lastIndexOf("/") + 1));
    FileUtils.copyResource(resource, samFile);
    if (indexed) {
      final String ext = resource.endsWith(".bam") ? ".bai" : ".tbi";
      FileUtils.copyResource(resource + ext, new File(samFile.getParentFile(), samFile.getName() + ext));
    }
    final Collection<File> coll = new ArrayList<>();
    coll.add(samFile);
    final SamReadingContext c = new SamReadingContext(coll, 1, filterParams, SamUtils.getUberHeader(coll), null);
    return new ThreadedMultifileIterator<>(c, new SingletonPopulatorFactory<>(new SamRecordPopulator()));
  }

  public static void check(SamFilterParams params, int expectedCount, int expectedFiltered, boolean indexed) throws IOException {
    try (TestDirectory dir = new TestDirectory("cnvproducttask")) {
      final RecordIterator<SAMRecord> multiIt = getIterator(dir, params, "com/rtg/variant/cnv/resources/testFilter" + (indexed ? ".bam" : ".sam"), indexed);
      try (SamFilterIterator it = new SamFilterIterator(multiIt, new DefaultSamFilter(params))) {
        int count = 0;
        while (it.hasNext()) {
          it.next();
          ++count;
        }
        assertEquals(expectedCount, count);
        assertEquals(expectedFiltered, it.getFilteredRecordsCount());
        assertEquals(false, it.hasNext());
        try {
          it.next();
          fail();
        } catch (final NoSuchElementException e) {
          //ignore
        }
      }
    }
  }

  static final String SAM_HEAD1 = ""
    + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n"
    + "@SQ" + TAB + "SN:gi" + TAB + "LN:30\n";
  static final String SAM_TAIL = TAB + "0" + TAB + "gi" + TAB + "%d" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS
    ;
  static final String SAM_TAIL2 = TAB + "0" + TAB + "gi2" + TAB + "%d" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + LS
    ;

  public void testFailureInRecords() throws IOException {
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      final StringBuilder sam1 = new StringBuilder();
      final StringBuilder sam2 = new StringBuilder();
      final StringBuilder sam3 = new StringBuilder();
      sam1.append(SAM_HEAD1);
      sam2.append(SAM_HEAD1);
      sam3.append(SAM_HEAD1);
      for (int i = 1; i < 21; ++i) {
        sam1.append("readA").append(i).append(String.format(SAM_TAIL, i));
        sam2.append("readB").append(i).append(String.format(SAM_TAIL, i));
        sam3.append("readC").append(i).append(String.format(SAM_TAIL, i));
        if (i == 10) {
          sam3.append("holy crap this file is broken").append(LS);
        }
      }
      final ArrayList<File> files = new ArrayList<>();
      int j = 0;
      for (final String content : new String[] {sam1.toString(), sam2.toString(), sam3.toString()}) {
        final File ff = new File(mDir, "alignments" + j + ".sam");
        FileUtils.stringToFile(content, ff);
        files.add(ff);
        ++j;
      }
      try {
        final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
        try (ThreadedMultifileIterator<SAMRecord> iterator = new ThreadedMultifileIterator<>(new SamReadingContext(files, 2, SamFilterParams.builder().create(), SamUtils.getUberHeader(files), null), pf)) {
          while (iterator.hasNext()) {
            final SAMRecord sr = iterator.next();
            //System.out.println(sr.getReadName());
            assertFalse(sr.getReadName().equals("readB12"));
          }
        }
        fail("Expected an exception");
      } catch (final NoTalkbackSlimException e) {
        //System.out.println(e.getMessage());
        assertEquals(SAMFormatException.class, e.getCause().getClass());
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  @SuppressWarnings("try")
  public void testFailureDuringSetup() throws IOException {
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      final ArrayList<File> files = new ArrayList<>();
      for (int j = 0; j < 10; ++j) {
        final StringBuilder sam = new StringBuilder(SAM_HEAD1);
        for (int i = 1; i < 21; ++i) {
          sam.append("readA").append(i).append(String.format(SAM_TAIL, i));
        }
        final File ffsam = new File(mDir, "alignments" + j + ".sam");
        FileUtils.stringToFile(sam.toString(), ffsam);
        files.add(ffsam);
      }
      final File incompatible = new File(mDir, "alignments-incompat.sam");
      FileUtils.stringToFile(SAM_HEAD1 + "readA-bad" + String.format(SAM_TAIL2, 22), incompatible);

      try {
        // This factory is just to induce some lag during the job scheduling (not the execution)
        final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<SAMRecord>(new SamRecordPopulator()) {
          @Override
          public Populator<SAMRecord> populator() {
            try {
              Thread.sleep(1);
            } catch (InterruptedException e) {
              // Ignore
            }
            return super.populator();
          }
        };
        final SAMFileHeader uberHeader = SamUtils.getUberHeader(files);
        assertTrue(SamUtils.checkHeaderDictionary(uberHeader, SamUtils.getSingleHeader(incompatible)));
        files.add(incompatible);
        try (ThreadedMultifileIterator<SAMRecord> ignored = new ThreadedMultifileIterator<>(new SamReadingContext(files, 2, SamFilterParams.builder().create(), uberHeader, null), pf)) {
          fail("This test needs to trigger an exception during iterator setup");
        }
        fail("This test needs to trigger an exception during iterator setup");
      } catch (final NoTalkbackSlimException e) {
        // expected exception about truncated sam record line
        assertTrue(e.getMessage().contains("Not enough fields")); // This must be the true cause of the failure, not a SlimAbortException that other workers produce when aborting.
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testZeroLengthFilesHandling() throws Exception {
    //Partitioning the files may end up with an empty file list - the two empty files will fill bin 2 and leave bin 3 untouched.
    final File f1 = new File(mDir, "one.sam");
    final File f2 = new File(mDir, "two.sam");
    final File f3 = new File(mDir, "three.sam");
    //first file needs a valid header.
    FileUtils.stringToFile(SAM_HEAD1, f1);
    assertTrue(f2.createNewFile());
    assertTrue(f3.createNewFile());
    final List<File> files = new ArrayList<>();
    files.add(f1);
    files.add(f2);
    files.add(f3);
    final SAMFileHeader sfh = new SAMFileHeader();
    sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    new ThreadedMultifileIterator<>(new SamReadingContext(files, 3, SamFilterParams.builder().create(), sfh, null), pf);

    //if we get here, it hasn't exploded. Success.
  }


  public void testFilterIterator() throws IOException {
    check(SamFilterParams.builder().create(), 22, 0, false);
    check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).create(), 12, 10, false);
    check(SamFilterParams.builder().maxAlignmentCount(0).maxMatedAlignmentScore(new IntegerOrPercentage(10)).create(), 2, 20, false);
    check(SamFilterParams.builder().maxAlignmentCount(1).maxMatedAlignmentScore(null).create(), 22, 0, false);
    check(SamFilterParams.builder().minMapQ(3).excludeUnmapped(true).create(), 17, 5, false);

    // These ones use a restriction which will remove some records prior to being seen by DefaultSamFilter

    check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).restriction("simulatedSequence1").create(), 6, 5, true); // Indexed results
    if (MultifileIterator.FALLBACK) { // Non-indexed uses fallback iterator
      check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).restriction("simulatedSequence1").create(), 6, 5, false);
    } else {
      try {
        check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).restriction("simulatedSequence1").create(), 6, 5, false);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("is not tabix indexed"));
      }

    }

    check(SamFilterParams.builder().maxAlignmentCount(-1).restriction("simulatedSequence1:90-145").create(), 6, 0, true);
  }

  public void testInvalidPopulate() throws IOException {

    try (TestDirectory dir = new TestDirectory()) {
      final File samFile = FileHelper.resourceToFile("com/rtg/sam/resources/badReadBase.sam", new File(dir, "samfile.sam"));

      final ThreadedMultifileIterator<VariantAlignmentRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(samFile), new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 2)));
      while (it.hasNext()) {
        it.next();
      }
      assertEquals(2L, it.getInvalidRecordsCount());
    }
  }

}
