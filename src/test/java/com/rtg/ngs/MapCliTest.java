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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class MapCliTest extends AbstractCliTest {

  private static final String APP_NAME = "rtg " + MapCli.MODULE_NAME;

  @Override
  protected AbstractCli getCli() {
    return new MapCli();
  }

  @Override
  public final void testApplicationName() {
    assertEquals(APP_NAME, new MapCli().applicationName() + " " + new MapCli().moduleName());
  }

  public final void testInitFlags() {
    checkHelp("maximum number of top equal results output per read (Default is 5",
        "output unmapped reads",
        "--" + MapFlags.NO_UNMATED, "do not output unmated reads when in paired-end mode",
        "--" + MapFlags.OUTPUT_UNFILTERED, "output all alignments meeting thresholds instead of applying mating and N limits"
        );
  }

  public void testFormatFlags() throws IOException {
    Diagnostic.setLogStream();
    final File mainOut = FileUtils.createTempDir("map", "test");
    try {
      final File left = new File(mainOut, "left");
      ReaderTestUtils.getReaderDNA(">l0\nacgt", left, null);
      final File right = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNA(">r0\nacgt", right, null);
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null);
      final File out = new File(mainOut, "out");
      checkHandleFlagsErr("");
      checkHandleFlagsErr("-i", mainOut.getPath());
      checkHandleFlagsErr("-i", mainOut.getPath(), "-t", template.getPath());
      checkHandleFlagsOut("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath());
      checkHandleFlagsErr("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath());
      checkHandleFlagsErr("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath());
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath());
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fasta"); //we don't validate file integrity in flag handling
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fastq");
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fastq", "-q", "bobo");
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fastq", "-q", "sanger");
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fastq", "-q", "solexa");
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-r", right.getPath(), "-F", "fastq", "-q", "illumina");
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-l", left.getPath(), "-F", "fastq", "-q", "sanger");
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-r", right.getPath(), "-F", "fastq", "-q", "sanger");
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-i", mainOut.getPath(), "-F", "fastq", "-q", "sanger");
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-i", mainOut.getPath(), "-F", "fasta", "-q", "sanger");
      checkHandleFlagsOut("-t", template.getPath(), "-o", out.getPath(), "-i", mainOut.getPath(), "-F", "fasta");
      checkHandleFlagsErr("-t", template.getPath(), "-o", out.getPath(), "-i", mainOut.getPath(), "-F", "bobo");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath(), "--sex", "male"), "is missing a 'reference.txt'");
    } finally {
      FileHelper.deleteAll(mainOut);
    }
  }

  public final void testMakeParamsCFlags() throws Exception {
    Diagnostic.setLogStream();
    final File mainOut = FileUtils.createTempDir("map", "test");
    try {
      final File left = new File(mainOut, "left");
      ReaderTestUtils.getReaderDNA(">l0\nacgt", left, null);
      final File right = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNA(">r0\nacgt", right, null);
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null);
      final File out = new File(mainOut, "out");
      final MapCli map = (MapCli) mCli;
      checkHandleFlagsOut("-i", mainOut.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath(),
          "--" + MapFlags.MIN_HITS_FLAG, "2");

      try (NgsParams params = map.makeParams()) {
        assertEquals(left.getPath(), params.buildFirstParams().directory().getPath());
        assertEquals(right.getPath(), params.buildSecondParams().directory().getPath());
        assertEquals(template.getPath(), params.searchParams().directory().getPath());
        assertEquals(out.getPath(), params.outputParams().directory().getPath());

        assertEquals(1000, params.maxFragmentLength().intValue());
        assertEquals(0, params.minFragmentLength().intValue());
        assertEquals(false, params.outputParams().progress());
        assertEquals(true, params.outputParams().isCompressOutput());
        assertEquals(false, params.outputParams().useids());
        assertEquals(false, params.outputParams().exclude());
        assertTrue(params.outputParams().outputUnmapped());

        assertEquals(false, params.useLongReadMapping());
        assertEquals(5, params.outputParams().topN());
        assertEquals(4095, params.outputParams().errorLimit());
        assertEquals(Integer.valueOf(2), params.minHits());
        assertEquals(-1, params.readFreqThreshold());

        //assertTrue(params.searchParams().reader() instanceof CompressedMemorySequencesReader2);

        assertEquals(out, map.outputDirectory());
        assertTrue(getCFlags().getUsageString().contains("Aligns sequence reads onto a reference template"));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(mainOut));
    }
  }


  public final void testGetTaskShortPE() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(bos);
    Diagnostic.setLogStream(ps);
    try {

      final MapCli map = (MapCli) mCli;

      final NgsParamsBuilder npb = new NgsParamsBuilder();
      final NgsOutputParams nop = new NgsOutputParams(new NgsOutputParamsBuilder());
      npb.outputParams(nop);
      final NgsParams params = new NgsParams(npb);

      final ParamsTask<?, ?> task = map.task(params, null);
      assertTrue(task instanceof NgsTask);
      ps.flush();

      assertEquals(0, task.usage());
      TestUtils.containsAll(bos.toString(), MapCli.MODULE_NAME + " paired=" + Boolean.FALSE + ", long=" + Boolean.FALSE + " running NgsTask");

    } finally {
      Diagnostic.setLogStream();
      ps.close();
    }
  }

  public void testFilterParams() {

    final CFlags flags = new CFlags();
    flags.registerOptional(MapFlags.TOPN_RESULTS_FLAG, Integer.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.MAX_TOP_RESULTS_FLAG, Integer.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.XSCORE_INDEL, Integer.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.REPORTING);

    flags.setFlags("--" + MapFlags.TOPN_RESULTS_FLAG, "5",
                          "--" + MapFlags.MAX_TOP_RESULTS_FLAG, "3",
                          "--" + MapFlags.XSCORE_INDEL, "4");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final NgsFilterParams filterparams = MapCli.makeFilterParams(flags, true, null);

      assertEquals(OutputFilter.TOPN_PAIRED_END, filterparams.outputFilter());

      assertTrue(filterparams.zip());
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
