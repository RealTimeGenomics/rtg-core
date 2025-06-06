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

package com.rtg.simulation;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;
import com.rtg.util.io.LogStream;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 */
public class ReadSimEvalCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ReadSimEvalCli();
  }

  private static final String READS = ReadSimEvalParamsTest.READS;

  private static final String SAM_ENDLINE = "\n";

  private static final String SAM = "@HD\tVN:1.0\tSO:coordinate" + SAM_ENDLINE + "@SQ\tSN:1\tLN:247249719" + SAM_ENDLINE
      + "read0:1:1234:S1:I0:D0\t16\t1\t1234\t255\t13=1X29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:1\tNM:i:0\tIH:i:1"
      + SAM_ENDLINE
      + "read1:1:2345R:S0:I1:D0\t16\t1\t23400\t255\t13=1X29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:2\tNM:i:0\tIH:i:1"
      + SAM_ENDLINE
     + "read2:1:56785:S4:I0:D0\t16\t1\t501889\t255\t12=2X29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:2\tNM:i:0\tIH:i:1"
      + SAM_ENDLINE;

  public final void testSingleEnd() throws IOException {
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READS, reads, null).close();
      final File sam = new File(topLevel, "samfile.sam.gz");
      FileHelper.stringToGzFile(SAM, sam);

      final MainResult res = checkMainInit("--output", rseDir.getPath(), "--reads", reads.getPath(), sam.getPath(), "--verbose", "--score-histogram");

      TestUtils.containsAll(res.err(),
          "No template map present in reads SDF. Cannot verify evaluation is against correct template." + LS,
          "Missing template SDF ID in SAM header. Cannot verify evaluation is against correct template." + LS
            );

      TestUtils.containsAll(res.out(),
          "Total SAM records = 3",
          "Total reads = 3",
          "Mapped reads = 3",
          "Unmapped reads = 0",
          "Reads mapped correctly = 1",
          "Reads mapped incorrectly = 2",
          "Accuracy (Precision) = correct / total mapped reads =  33.33", "Sensitivity (Recall) = correct / total reads =  33.33",
          "Reads mapped correctly and unique = 1",
          "Reads mapped correctly and multiple = 0",
          "Reads mapped incorrectly and unique = 2",
          "Reads mapped incorrectly and multiple = 0"
          );

      TestUtils.containsAll(FileUtils.fileToString(new File(rseDir, "score_hist.tsv")),
          "Alignment Score Distribution",
          "Score\ttrue_positives\tfalse_positives",
          //"0\t0\t0",
          "1\t1.00\t0.00",
          "2\t0.00\t2.00",
          "Generated Score Distribution",
          "Score\ttrue_positives\tfalse_positives",
          //"0\t0\t0",
          "1\t1.00\t0.00",
          "2\t0.00\t1.00",
          //"3\t0\t0",
          "4\t0.00\t1.00"
          );
    }
  }

  private static final String SAM_PAIRED = "@HD\tVN:1.0\tSO:coordinate" + SAM_ENDLINE
      + "@SQ\tSN:1\tLN:247249719" + SAM_ENDLINE
      //                                                             0123456789012345678901234567890123456789012
      + "read0:1:1234:S1:I0:D0\t89\t1\t1234\t255\t13=1X29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:1\tNM:i:0\tIH:i:1" + SAM_ENDLINE
      + "read1:1:2345R:S0:I1:D0\t89\t1\t23400\t255\t13=1X29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:2\tNM:i:0\tIH:i:1" + SAM_ENDLINE
      + "read2:1:56785:S4:I0:D0\t89\t1\t501889\t255\t14=2N29=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:2\tNM:i:0\tIH:i:1" + SAM_ENDLINE
      ;

  /**
   * Test method for
   * {@link ReadSimEvalCli#mainExec(java.io.OutputStream, java.io.PrintStream)}
   * .
   */
  public final void testPairedEnd() throws IOException {
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READS, new File(reads, "left"), null).close();
      ReaderTestUtils.getReaderDNA(READS, new File(reads, "right"), null).close();
      final File sam = new File(topLevel, "samfile.sam.gz");
      FileHelper.stringToGzFile(SAM_PAIRED, sam);

      final MainResult res = checkMainInit("--output", rseDir.getPath(), "--reads", reads.getPath(), sam.getPath(), "--verbose", "--score-histogram", "--mapq-histogram");

      TestUtils.containsAll(res.out(),
            "Total SAM records = 3",
            "Total pairs = 3",
            "Mated pairs = 0",
            "Unmated pairs = 3",
            "Left reads mapped correctly = 1",
            "Left reads mapped incorrectly = 2",
            "Right reads mapped correctly = 0",
            "Right reads mapped incorrectly = 0",
            "Accuracy (Precision) = correct / total mapped reads =  33.33",
            "Sensitivity (Recall) = correct / total reads =  16.67",
            "Left reads mapped correctly and unique = 1",
            "Left reads mapped correctly and multiple = 0",
            "Left reads mapped incorrectly and unique = 2",
            "Left reads mapped incorrectly and multiple = 0",
            "Right reads mapped correctly and unique = 0",
            "Right reads mapped correctly and multiple = 0",
            "Right reads mapped incorrectly and unique = 0",
            "Right reads mapped incorrectly and multiple = 0"
            );
      TestUtils.containsAll(FileUtils.fileToString(new File(rseDir, "mapq_hist.tsv")),
            "MAPQ Score Distribution",
            "Score\ttrue_positives\tfalse_positives",
            "255\t1.00\t2.00"
            );
      TestUtils.containsAll(FileUtils.fileToString(new File(rseDir, "score_hist.tsv")),
        "Alignment Score Distribution",
        "Score\ttrue_positives\tfalse_positives",
        //"0\t0\t0",
        "1\t1.00\t0.00",
        "2\t0.00\t2.00",
        "Generated Score Distribution",
        "Score\ttrue_positives\tfalse_positives",
        //"0\t0\t0",
        "1\t1.00\t0.00",
        "2\t0.00\t1.00",
        //"3\t0\t0",
        "4\t0.00\t1.00"
      );
    }
  }

  private static final String READ_LEFT = ">read0/0/0/simulatedSequence1/27/R/10." + LS
                                        + "TTATACACGA" + LS
                                        + ">read1/0/0/simulatedSequence1/69/F/10." + LS
                                        + "TATTAAGCTC" + LS
                                        + ">read2/0/0/simulatedSequence1/66/F/10." + LS
                                        + "GACTATTAAG" + LS;

  private static final String READ_RIGHT = ">read0/0/0/simulatedSequence1/14/F/10." + LS
                                        + "TCCCAACTTT" + LS
                                        + ">read1/0/0/simulatedSequence1/80/R/10." + LS
                                        + "TGGGCCGTCT" + LS
                                        + ">read2/0/0/simulatedSequence1/80/R/10." + LS
                                        + "TGGGCCGTCT" + LS;


  private static final String SAM_PAIRED2 = "@HD VN:1.0 SO:coordinate" + SAM_ENDLINE
      + "@SQ SN:simulatedSequence1 LN:100" + SAM_ENDLINE
      + "read0/0/0/simulatedSequence1/14/F/10. 163 simulatedSequence1 14 55 10= = 27 23 TCCCAACTTT 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE
      + "read0/0/0/simulatedSequence1/27/R/10. 83 simulatedSequence1 27 55 10= = 14 -23 TCGTGTATAA 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE
      + "read2/0/0/simulatedSequence1/66/F/10. 99 simulatedSequence1 66 55 10= = 80 24 GACTATTAAG 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE
      + "read1/0/0/simulatedSequence1/69/F/10. 99 simulatedSequence1 69 55 10= = 80 21 TATTAAGCTC 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE
      + "read1/0/0/simulatedSequence1/80/R/10. 147 simulatedSequence1 80 55 10= = 69 -21 AGACGGCCCA 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE
      + "read2/0/0/simulatedSequence1/80/R/10. 147 simulatedSequence1 80 55 10= = 66 -24 AGACGGCCCA 5555555555 AS:i:0 NM:i:0 MQ:i:255 XA:i:0 IH:i:1 NH:i:1" + SAM_ENDLINE;

  public final void testPairedEndMated() throws IOException {
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READ_LEFT, new File(reads, "left"), null).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, new File(reads, "right"), null).close();
      final File sam = new File(topLevel, "samfile.sam.gz");
      FileHelper.stringToGzFile(SAM_PAIRED2.replaceAll(" ", "\t"), sam);

      final MainResult res = checkMainInit("--output", rseDir.getPath(), "--reads", reads.getPath(), sam.getPath(), "--verbose");
      //System.err.println(bos.toString());
      TestUtils.containsAll(res.out(),
        "Total SAM records = 6",
        "Total pairs = 3",
        "Mated pairs = 3",
        "Unmated pairs = 0",
        "Left reads mapped correctly = 3",
        "Left reads mapped incorrectly = 0",
        "Right reads mapped correctly = 3",
        "Right reads mapped incorrectly = 0",
        "Accuracy (Precision) = correct / total mapped reads = 100.00",
        "Sensitivity (Recall) = correct / total reads = 100.00"
      );
    }
  }

  private static final String SAM_NO_NAME = "@HD\tVN:1.0\tSO:coordinate" + SAM_ENDLINE + "@SQ\tSN:1\tLN:247249719" + SAM_ENDLINE
  //                                             0123456789012345678901234567890123456789012
      + "read0\t89\t1\t1234\t255\t22=1X20=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:1\tNM:i:0\tIH:i:1" + SAM_ENDLINE;

  public void testErr() throws IOException {
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READS, new File(reads, "left"), null).close();
      ReaderTestUtils.getReaderDNA(READS, new File(reads, "right"), null).close();
      final File sam = new File(topLevel, "samfile.sam.gz");
      FileHelper.stringToGzFile(SAM_NO_NAME, sam);

      final String err = checkMainInitBadFlags("--output", rseDir.getPath(), "--reads", reads.getPath(), sam.getPath());
      assertTrue(err.contains("does not appear to be"));
    }
  }

  private static final String SAM_BROKEN_READNAME = "@HD\tVN:1.0\tSO:coordinate" + SAM_ENDLINE
  + "@SQ\tSN:1\tLN:247249719" + SAM_ENDLINE
  + "stupid_read_name\t16\t1\t1234\t255\t22=1X20=\t*\t0\t0\tACATGCTGCATGCATGCTGATGCTGCTGCTATAGTGATGAATC\t*\tAS:i:1\tNM:i:0\tIH:i:1"
  + SAM_ENDLINE;

  public void testWarning() throws IOException {
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      final File rseDir = new File(topLevel, "rse");
      ReaderTestUtils.getReaderDNA(READS, reads, null).close();
      final File sam = new File(topLevel, "samfile.sam.gz");
      FileHelper.stringToGzFile(SAM_BROKEN_READNAME, sam);

      final String err = checkMainInitBadFlags("--output", rseDir.getPath(), "--reads", reads.getPath(), sam.getPath(), "--verbose");
      assertTrue(err.contains("Read stupid_read_name does not appear to be a simulated read"));
    }
  }

  public void testSoftClippingAtStart() throws Exception {
    final ReadSimEvalCli rma = new ReadSimEvalCli() {
      @Override
      protected int mainExec(OutputStream out, LogStream log) {
        return 0;
      }
    };
    try (final TestDirectory topLevel = new TestDirectory("readmappingaccuracy")) {
      final File reads = new File(topLevel, "reads");
      ReaderTestUtils.getReaderDNA(">read0/0/1/chr1/4730076/R/4.1D8.1D55.1D27.1D6." + LS
          + "TTCCATGGTTTTCTGAGCCTCAGTTTTCTCATCTGAAAAACGGGGATGTCACTCAGCCCTGCACAGGCTGGAAGGATGGTGACCCCCTACCATTACAGGT" + LS, reads, null).close();

      final MemoryPrintStream recordsMps = new MemoryPrintStream();
      try (LineWriter recordslw = new LineWriter(new OutputStreamWriter(recordsMps.outputStream()))) {
        final MemoryPrintStream mps = new MemoryPrintStream();
        final File f = FileHelper.createTempFile(topLevel);
        rma.mainInit(new String[]{"--reads", reads.getPath(), "-o", new File(topLevel, "out").getPath(), f.getPath()}, mps.outputStream(), mps.printStream());
        rma.init();

        final SAMFileHeader sfh = new SAMFileHeader();

        //read0/0/1/chr1/4730076/R/4.1D8.1D55.1D27.1D6.        16      chr1    4730086 71      8S57M1D27M1D8M  *       0       0       TTCCATGGTTTTCTGAGCCTCAGTTTTCTCATCTGAAAAACGGGGATGTCACTCAGCCCTGCACAGGCTGGAAGGATGGTGACCCCCTACCATTACAGGT    5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555    RG:Z:NOID       PG:Z:tmap       MD:Z:57^G27^T8  NM:i:2  AS:i:78 XA:Z:map3-1     XS:i:-2147483647        XT:i:10 ZS:i:4729967    ZE:i:4730267
        final SAMRecord rec = new SAMRecord(sfh);
        rec.setAlignmentStart(4730086);
        rec.setReadName("read0/0/1/chr1/4730076/R/4.1D8.1D55.1D27.1D6.");
        rec.setCigarString("8S57M1D27M1D8M");
        rec.setFlags(16);
        rec.setReferenceName("chr1");

        rma.processRecord(rec, recordslw);

        assertTrue(rma.mLeftStats.isMapped(0));
        assertFalse(rma.mLeftStats.isMultiple(0));
        assertTrue(rma.mLeftStats.isFound(0));
      }
    }
  }

}
