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
package com.rtg.ngs;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.tempstage.PairedTempFileWriterImpl;
import com.rtg.pairedend.SlidingWindowCollector;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class UnfilteredPairedEndOutputProcessorTest extends AbstractPairedEndOutputProcessorSyncTest {

  @Override
  OutputProcessor getPairedEndOutputProcessorSync(NgsParams param, MapStatistics stats, boolean outputUnmated, boolean outputUnmapped)
      throws IOException {
    return new UnfilteredPairedEndOutputProcessor(param, stats, outputUnmapped);
  }

  @Override
  PairedEndOutputProcessor getPEOP(PairedTempFileWriterImpl writer, SlidingWindowCollector collector) {
    return new PairedEndOutputProcessor(writer, collector);
  }

  @Override
  String getThreadCloneExpectedResource() {
    return "tnpeops_pe_unfiltered.txt";
  }

  @Override
  OutputFilter getOutputFilter() {
    return OutputFilter.SAM_UNFILTERED;
  }

  @Override
  String getOutputFilePrefix() {
    return NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME;
  }
  @Override
  String getOutputBamFilePrefix() {
    return NgsOutputParams.ALIGNMENTS_BAM_FILE_NAME;
  }

  public void testUnmatedLocal() throws Exception {
    final File tmp = FileUtils.createTempDir("topnsync", "clone");
    try {
      checkUnmated(tmp);
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testUnmatedRemote() throws Exception {
    checkUnmated(null);
  }

  /**
   * Test of threadClone method, of class PairedEndOutputProcessorSync.
   */
  public void checkUnmated(final File tempDir) throws Exception {
    final ByteArrayOutputStream logBytes = new ByteArrayOutputStream();
    try (PrintStream log = new PrintStream(logBytes)) {
      Diagnostic.setLogStream(log);
      final int numThreads = 4;
      final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "TestUnmated", true);
      final NgsParams params = getDefaultBuilder(tempDir, false, OutputFilter.SAM_UNFILTERED, null).numberThreads(numThreads).create();
      try (UnfilteredPairedEndOutputProcessor sync = new UnfilteredPairedEndOutputProcessor(params, null, true)) {
        for (int i = 0; i < numThreads; ++i) {
          final long padding = params.calculateThreadPadding();
          final long start = i * MAX_COORD / numThreads;
          final long end = (i + 1) * MAX_COORD / numThreads;
          final HashingRegion region = new HashingRegion(0, start, 0, end, Math.max(0, start - padding), Math.min(MAX_COORD, end + padding));
          stp.execute(new SimpleProcess2(sync, region, i));
        }
        stp.terminate();
        sync.finish();
      }
      final String contentsUnmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMAPPED_SAM_FILE_NAME)));
      final String contentsUnmated = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME)));
      mNano.check("upeops-unmated", contentsUnmated, false);
      mNano.check("upeops-unmapped", contentsUnmapped, false);
    }
  }

  public void testSortedAppend() throws IOException {
    sortedAppendTest("tnpeops_pe_unfiltered_sortedAppend.txt");
  }
  public void testSortedAppendBam() throws IOException {
    sortedAppendTestBam("tnpeops_pe_unfiltered_sortedAppend.txt");
  }
  //                                  0    0    1    1    2    2    3    3
  //                                  0    5    0    5    0    5    0    5/
  static final String TEMPLATE_STR = "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg";
  static final String TEMPLATE = ">t" + StringUtils.LS + TEMPLATE_STR + StringUtils.LS;
  static final String LEFT_READ_AS2 =  "acacactgcaagcaagagggcctccc";            //starts at 0
  static final String LEFT_READ_AS6 =  "acacactgcggtctgagagggcctcccac";
  static final String RIGHT_READ_AS0 = "cccctttggcccccgaccagtgtgggctga";        //starts at 36
  static final String RIGHT_READ_AS8 = "cccctttttaaaaagagcagtgtgggctga";
  static final String LEFT_READ_AS11 = "ctgcaactgttctaaagctcccacagcactct";      //starts at tpos 5
  static final String RIGHT_READ_AS11 = "agagtgctgtgggagctttagaacagttgcag";
  static final String READ_LEFT = ">r0" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r1" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r2" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r3" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r4" + StringUtils.LS + LEFT_READ_AS6 + StringUtils.LS
                          + ">r5" + StringUtils.LS + LEFT_READ_AS11 + StringUtils.LS
                          + ">r6" + StringUtils.LS + LEFT_READ_AS2 + StringUtils.LS
                          + ">r7" + StringUtils.LS + LEFT_READ_AS11 + StringUtils.LS;
  static final String READ_RIGHT = ">r0" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS
                          + ">r1" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS
                          + ">r2" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r3" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r4" + StringUtils.LS + RIGHT_READ_AS11 + StringUtils.LS
                          + ">r5" + StringUtils.LS + RIGHT_READ_AS8 + StringUtils.LS
                          + ">r6" + StringUtils.LS + RIGHT_READ_AS11 + StringUtils.LS
                          + ">r7" + StringUtils.LS + RIGHT_READ_AS0 + StringUtils.LS;

  public void testEndToEnd() throws Exception {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final File template = FileUtils.createTempDir("template", "ngs", tmpDir);
      final File left = new File(tmpDir, "left");
      final File right = new File(tmpDir, "right");
      final File out = new File(tmpDir, "out");

      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

      try (MemoryPrintStream mps = new MemoryPrintStream()) {
        final MapCli map = new MapCli();
        assertEquals(0, map.mainInit(new String[]{
          "-e", "10",
          "-E", "7",
          "--all-hits",
          "-w", "4",
          "-s", "1",
          "-m", "1",
          "-i", tmpDir.toString(),
          "-t", template.toString(),
          "-o", out.toString(),
          "--no-gzip",
          "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "1",
          "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "1",
          "--" + MapFlags.MISMATCH_PENALTY_FLAG, "1",
          "--" + MapFlags.ALIGNER_MODE_FLAG, "general",
          "--" + MapFlags.SAM_FLAG,
          "--" + MapFlags.DONT_UNIFY_FLAG,
          "--" + MapFlags.UNKNOWNS_PENALTY_FLAG, "1"}, mps.outputStream(), mps.printStream()));

        final String alignments = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(out, "alignments.sam")));
        final String unmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(out, "unmapped.sam")));
        mNano.check("upeops-endtoend-alignments", alignments, false);
        mNano.check("upeops-endtoend-unmapped", unmapped, false);
      }
    }
  }

  public void testEndToEnd2() throws Exception {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final File template = FileUtils.createTempDir("template", "ngs", tmpDir);
      final File left = new File(tmpDir, "left");
      final File right = new File(tmpDir, "right");
      final File out = new File(tmpDir, "out");

      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

      final MemoryPrintStream mps = new MemoryPrintStream();

      final MapFCli map = new MapFCli();
      final int code = map.mainInit(new String[] {"-e", "10", "-E", "7", "-w", "4", "-s", "1", "-m", "1", "-i", tmpDir.toString(), "-t", template.toString(), "-o", out.toString()}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);

      assertFalse(new File(out, "alignments.sam").isFile());
      assertFalse(new File(out, "unmapped.sam").isFile());
      assertTrue(new File(out, "alignments.sdf").isDirectory());
      assertTrue(new File(out, "unmapped.sdf").isDirectory());
    }
  }

}
