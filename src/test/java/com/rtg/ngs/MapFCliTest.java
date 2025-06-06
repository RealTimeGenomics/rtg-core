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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.ngs.MapFCli.MapFTask;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.reader.SdfId;
import com.rtg.usage.UsageMetric;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 *
 */
public class MapFCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MapFCli();
  }

  public final void testInitFlags() {
    checkHelp("directory for output",
        "guaranteed number of positions",
        "guaranteed minimum number of indels",
        "reads",  // for the single-end -i flag.
        "maximum mismatches for mappings in single-end mode",
        "maximum mismatches for mappings across mated results",
        "maximum mismatches for mappings of unmated results",
        "do not gzip the output",
        "maximum repeat frequency",
        "guaranteed minimum number of substitutions",
        "number of threads",
        "number of available cores",
        "word size",
        "template",
        "Filters reads for contaminant sequences by mapping them against the contaminant reference.",
        "It outputs two SDF files, one containing the input reads that map to the reference and one containing those that do not."
        );
  }

  public final void testMapF() throws InvalidParamsException, IOException {

    final File tempDir = FileUtils.createTempDir("temp", "dir");
    try {
      final MapFCli mf = new MapFCli();
      assertEquals("mapf", mf.moduleName());
      final File outputFile = new File(tempDir, "testoutput");
      final File inputFile = new File(tempDir, "input");
      ReaderTestUtils.getReaderDNA(">x\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt", inputFile, new SdfId(5L));

      MainResult.run(mf, "-o", outputFile.getPath(), "-i", inputFile.getPath(), "-t", inputFile.getPath(), "--Xforce-long");
      assertEquals("testoutput", mf.outputDirectory().getName());

      final NgsParams params = mf.makeParams();
      new DummyTask(params, TestUtils.getNullOutputStream()).assertInstanceSE();
      assertEquals(11, params.stepSize());
      assertTrue(params.compressHashes());
      assertTrue(params.outputParams().filter().zip());
      assertEquals(1, params.readFreqThreshold());

      FileHelper.deleteAll(outputFile);

      //Now can have step > window
      checkMainInitOk("-o", outputFile.getPath(), "-i", inputFile.getPath(), "-t", inputFile.getPath(), "--Xforce-long", "--step", "23");
      final String usageLog = mf.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=mapf runId=", ", Usage end module=mapf runId=", " metric=116 success=true]");
    } finally {
      FileHelper.deleteAll(tempDir);
    }
  }

  public final void testRGMapF() throws InvalidParamsException, IOException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final MapFCli mf = new MapFCli();
      assertEquals("mapf", mf.moduleName());
      final File outputFile = new File(tempDir, "testoutput");
      final File inputFile = new File(tempDir, "input");
      ReaderTestUtils.getReaderDNA(">x\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt", inputFile, new SdfId(5L));
      MainResult.run(mf,
        "-o", outputFile.getPath(),
        "-i", inputFile.getPath(), "-t", inputFile.getPath(), "--sam-rg", "@RG\\tPL:IONTORRENT\\tSM:foo\\tID:foo");
      assertEquals("testoutput", mf.outputDirectory().getName());
      final NgsParams params = mf.makeParams();
      new DummyTask(params, TestUtils.getNullOutputStream()).assertInstanceSE();
      assertEquals("IONTORRENT", params.outputParams().readGroup().getPlatform());
      final String usageLog = mf.usageLog();
//      System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=mapf runId=", ", Usage end module=mapf runId=", " metric=116 success=true]");
    }
  }

  private static class DummyTask extends MapFTask {
    DummyTask(NgsParams params, OutputStream defaultOutput) {
      super(params, defaultOutput, new UsageMetric());
    }
    public void assertInstanceSE() {
      assertTrue(mStatistics instanceof MapFilterSingleEndMapStatistics);
    }
  }

  static final String TEMPLATE_STR = "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg";
  static final String TEMPLATE = ">t\n" + TEMPLATE_STR + "\n";
  static final String LEFT_READ_AS2 =  "ACACACTGCAAGCAAGAGGGCCTCCC";            //starts at 0
  static final String RIGHT_READ_AS0 = "CCCCTTTGGCCCCCGACCAGTGTGGGCTGA";        //starts at 36
  static final String LEFT_READ_AS6 =  "ACACACTGCGGTCTGAGAGGGCCTCCCAC";
  static final String RIGHT_READ_AS8 = "CCCCTTTTTAAAAAGAGCAGTGTGGGCTGA";
  static final String LEFT_READ_AS11 = "CTGCAACTGTTCTAAAGCTCCCACAGCACTCT";      //starts at tpos 5
  static final String RIGHT_READ_AS11 = "AGAGTGCTGTGGGAGCTTTAGAACAGTTGCAG";
  static final String READ_LEFT = ">r0\n" + LEFT_READ_AS2 + "\n"
      + ">r1\n" + LEFT_READ_AS6 + "\n"
      + ">r2 a\n" + LEFT_READ_AS6 + "\n"
      + ">r3\n" + LEFT_READ_AS2 + "\n"
      + ">r4 b\n" + LEFT_READ_AS6 + "\n"
      + ">r5\n" + LEFT_READ_AS11 + "\n"
      + ">r6 c\n" + LEFT_READ_AS2 + "\n"
      + ">r7\n" + LEFT_READ_AS11 + "\n";
  static final String READ_RIGHT = ">r0\n" + RIGHT_READ_AS0 + "\n"
      + ">r1\n" + RIGHT_READ_AS0 + "\n"
      + ">r2 a\n" + RIGHT_READ_AS8 + "\n"
      + ">r3\n" + RIGHT_READ_AS8 + "\n"
      + ">r4 b\n" + RIGHT_READ_AS11 + "\n"
      + ">r5\n" + RIGHT_READ_AS8 + "\n"
      + ">r6 c\n" + RIGHT_READ_AS11 + "\n"
      + ">r7\n" + RIGHT_READ_AS0 + "\n";


  public void testEndToEnd() throws Exception {
    Diagnostic.setLogStream();

    try (final TestDirectory tmpDir = new TestDirectory("mapftest")) {
      final File left = new File(tmpDir, "left");
      final File right = new File(tmpDir, "right");
      final File template = FileUtils.createTempDir("template", "ngs", tmpDir);
      final File out = new File(tmpDir, "out");

      ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();
      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();

      checkMainInitOk("--sam",
        "-e", "10",
        "-E", "7",
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
        "--" + MapFlags.DONT_UNIFY_FLAG,
        "--" + MapFlags.UNKNOWNS_PENALTY_FLAG, "1",
        "--" + MapFlags.ALIGNER_MODE_FLAG, "general"
      );

      final String alignments = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(out, "alignments.sam")));
      final String unmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(out, "unmapped.sam")));
      mNano.check("mapf-endtoend-alignments", alignments, false);
      mNano.check("mapf-endtoend-unmapped", unmapped, false);
      MainResult res = MainResult.run(new Sdf2Fasta(), "-i", new File(out, "alignments.sdf").toString(), "-o", new File(out, "outfa").toString(), "-Z");
      assertEquals(res.err(), 0, res.rc());

      assertEquals(">r0\n" + LEFT_READ_AS2 + "\n"
        + ">r1\n" + LEFT_READ_AS6 + "\n"
        + ">r2 a\n" + LEFT_READ_AS6 + "\n"
        + ">r3\n" + LEFT_READ_AS2 + "\n"
        + ">r4 b\n" + LEFT_READ_AS6 + "\n"
        + ">r6 c\n" + LEFT_READ_AS2 + "\n"
        + ">r7\n" + LEFT_READ_AS11 + "\n", FileUtils.fileToString(new File(out, "outfa_1.fasta")));
      assertEquals(">r0\n" + RIGHT_READ_AS0 + "\n"
        + ">r1\n" + RIGHT_READ_AS0 + "\n"
        + ">r2 a\n" + RIGHT_READ_AS8 + "\n"
        + ">r3\n" + RIGHT_READ_AS8 + "\n"
        + ">r4 b\n" + RIGHT_READ_AS11 + "\n"
        + ">r6 c\n" + RIGHT_READ_AS11 + "\n"
        + ">r7\n" + RIGHT_READ_AS0 + "\n", FileUtils.fileToString(new File(out, "outfa_2.fasta")));

      res = MainResult.run(new Sdf2Fasta(), "-i", new File(out, "unmapped.sdf").toString(), "-o", new File(out, "unmappedfa").toString(), "-Z");
      assertEquals(res.err(), 0, res.rc());

      assertEquals(">r5\n" + LEFT_READ_AS11 + "\n", FileUtils.fileToString(new File(out, "unmappedfa_1.fasta")));
      assertEquals(">r5\n" + RIGHT_READ_AS8 + "\n", FileUtils.fileToString(new File(out, "unmappedfa_2.fasta")));
    }
  }
}
