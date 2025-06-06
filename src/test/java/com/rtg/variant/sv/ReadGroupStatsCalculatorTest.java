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
package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.AbstractTest;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class ReadGroupStatsCalculatorTest extends AbstractTest {

  /** a SAM file string */
  public static final String RG_SAM = ""
    + "@HD\tVN:1.0\tSO:coordinate" + StringUtils.LS
    + "@RG\tID:boo_pe\tPL:ILLUMINA\tSM:bar" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence1\tLN:100000" + StringUtils.LS
    + "144\t163\tsimulatedSequence1\t27\t255\t35=\t=\t204\t212\tCTAAGCACTCAAGCTGGAGATTACCATACTTAGGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1371\t99\tsimulatedSequence1\t76\t255\t35=\t=\t254\t213\tGCGCTCGTAAATTCTCGACATTCCGCAGTGGCAGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1283\t99\tsimulatedSequence1\t130\t255\t35=\t=\t272\t177\tTCGCTGTGGAACTAGACCCGCCCTACTGTTGTCGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "2502\t163\tsimulatedSequence1\t131\t255\t35=\t=\t331\t235\tCGCTGTGGAACTAGACCCGCCCTACTGTTGTCGCT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "144\t83\tsimulatedSequence1\t204\t255\t35=\t=\t27\t-212\tTTTCGAATAGCTGGCGTGCACGATTGGCGTGTACA\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "1371\t147\tsimulatedSequence1\t254\t255\t35=\t=\t76\t-213\tACAGTTTCGTAGGGTTTCGTAGCAATGGCTTTGTC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "1061\t163\tsimulatedSequence1\t264\t255\t35=\t=\t457\t228\tAGGGTTTCGTAGCAATGGCTTTGTCTGAAGAGAAG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "1283\t147\tsimulatedSequence1\t272\t255\t35=\t=\t130\t-177\tGTAGCAATGGCTTTGTCTGAAGAGAAGGGGTCTGT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1890\t99\tsimulatedSequence1\t322\t255\t35=\t=\t426\t139\tCCCCAAAGTACCAGCGGGTGATCTCCGCAGACCTC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "2502\t83\tsimulatedSequence1\t331\t255\t35=\t=\t131\t-235\tACCAGCGGGTGATCTCCGCAGACCTCGGAGCCGAG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "218\t163\tsimulatedSequence1\t390\t255\t35=\t=\t539\t184\tTTAGGTCGGAGCAAGCGGCGTCCGTTAGAAAGCGT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1890\t147\tsimulatedSequence1\t426\t255\t35=\t=\t322\t-139\tGCTACTAGCCCGCGAGAGGCGACGTCTCACTCCCC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "640\t163\tsimulatedSequence1\t447\t255\t35=\t=\t633\t221\tACGTCTCACTCCCCGTAACAGCAGGACTGGTCTGG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1061\t83\tsimulatedSequence1\t457\t255\t35=\t=\t264\t-228\tCCCCGTAACAGCAGGACTGGTCTGGGTACGGCTGG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1790\t99\tsimulatedSequence1\t458\t255\t35=\t=\t540\t117\tCCCGTAACAGCAGGACTGGTCTGGGTACGGCTGGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "99\t99\tsimulatedSequence1\t461\t255\t35=\t=\t569\t143\tGTAACAGCAGGACTGGTCTGGGTACGGCTGGCGCA\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1048\t99\tsimulatedSequence1\t483\t255\t35=\t=\t708\t260\tTACGGCTGGCGCAGGCATGAATGGGATACCCGTGT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "57\t99\tsimulatedSequence1\t517\t255\t35=\t=\t726\t244\tTTCTGCCTTCTGTGTCTACATCCGGTTATGGGAGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1243\t163\tsimulatedSequence1\t539\t255\t35=\t=\t720\t216\tCGGTTATGGGAGCTAGCGAAAGTCACGCGAACCAG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "218\t83\tsimulatedSequence1\t539\t255\t35=\t=\t390\t-184\tCGGTTATGGGAGCTAGCGAAAGTCACGCGAACCAG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1790\t147\tsimulatedSequence1\t540\t255\t35=\t=\t458\t-117\tGGTTATGGGAGCTAGCGAAAGTCACGCGAACCAGT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1867\t99\tsimulatedSequence1\t562\t255\t35=\t=\t688\t161\tCACGCGAACCAGTAAGGTTAGAGTAGAAACGGAGT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "99\t147\tsimulatedSequence1\t569\t255\t35=\t=\t461\t-143\tACCAGTAAGGTTAGAGTAGAAACGGAGTACTACTT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "1178\t99\tsimulatedSequence1\t606\t255\t35=\t=\t743\t172\tCTATTGAAAAAAGGTGCTGCATCTGAGCCACGCCT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    + "640\t65\tsimulatedSequence1\t633\t255\t35=\t=\t447\t-221\tCCACGCCTAATTGCTGAGCTCTAAGGCTGTGAGTG\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tRG:Z:boo_pe\tIH:i:1" + StringUtils.LS
    ;

  public void testReadGroupStatsCalculator() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory("readgroup")) {
      final File samFile = FileUtils.stringToFile(RG_SAM, new File(tempDir, "rg.sam"));
      final List<File> files = new ArrayList<>();
      files.add(samFile);
      final ReadGroupStatsCalculator calculator = new ReadGroupStatsCalculator();
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      calculator.calculate(files, out.outputStream());
      //System.err.println(out.toString());
      //System.err.println(err.toString());
      TestUtils.containsAll(out.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "boo_pe\t"
          );
      TestUtils.containsAll(err.toString(),
          "Accumulating read group statistics",
          "Skipping record with no RG id",
          "Skipped 6 records with no RG id"
          );
    }
  }

  static final String NO_RG_SAM = ""
    + "@HD\tVN:1.0\tSO:coordinate" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence1\tLN:400" + StringUtils.LS
    + "144\t163\tsimulatedSequence1\t27\t255\t35=\t=\t204\t212\tCTAAGCACTCAAGCTGGAGATTACCATACTTAGGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "1371\t99\tsimulatedSequence1\t76\t255\t35=\t=\t254\t213\tGCGCTCGTAAATTCTCGACATTCCGCAGTGGCAGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "1283\t99\tsimulatedSequence1\t130\t255\t35=\t=\t272\t177\tTCGCTGTGGAACTAGACCCGCCCTACTGTTGTCGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "2502\t163\tsimulatedSequence1\t131\t255\t35=\t=\t331\t235\tCGCTGTGGAACTAGACCCGCCCTACTGTTGTCGCT\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    + "144\t83\tsimulatedSequence1\t204\t255\t35=\t=\t27\t-212\tTTTCGAATAGCTGGCGTGCACGATTGGCGTGTACA\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    ;

  public void testReadGroupStatsCalculatorNoReadGroups() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory("readgroup")) {
      final File samFile = FileHelper.stringToGzFile(NO_RG_SAM, new File(tempDir, "rg.sam.gz"));
      final List<File> files = new ArrayList<>();
      files.add(samFile);
      CommandLine.clearCommandArgs();
      final ReadGroupStatsCalculator calculator = new ReadGroupStatsCalculator();
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      try {
        calculator.calculate(files, out.outputStream());
      } catch (NoTalkbackSlimException e) {
        // Expected
      }
    }
  }

  static final String UNPAIRED_SAM = ""
    + "@HD\tVN:1.0\tSO:coordinate" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence1\tLN:400" + StringUtils.LS
    + "144\t0\tsimulatedSequence1\t27\t255\t35=\t=\t204\t212\tCTAAGCACTCAAGCTGGAGATTACCATACTTAGGC\t+++++++++++++++++++++++++++++++++++\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1" + StringUtils.LS
    ;

  public void testNonPairedSam() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory("readgroup")) {
      final File samFile = FileHelper.stringToGzFile(UNPAIRED_SAM, new File(tempDir, "rg.sam.gz"));
      final List<File> files = new ArrayList<>();
      files.add(samFile);
      final ReadGroupStatsCalculator calculator = new ReadGroupStatsCalculator();
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      calculator.calculate(files, NullStreamUtils.getNullOutputStream());
      TestUtils.containsAll(err.toString(), "No properly mated alignments were available");
    }
  }

  public void testMerge() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory("readgroup")) {
      final File samFile = FileUtils.stringToFile(RG_SAM, new File(tempDir, "rg.sam"));
      final List<File> files = new ArrayList<>();
      files.add(samFile);
      final ReadGroupStatsCalculator calculator = new ReadGroupStatsCalculator();
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      calculator.calculate(files, out.outputStream());
      //System.err.println(out.toString());
      //System.err.println(err.toString());
      TestUtils.containsAll(out.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "boo_pe\t"
          );
      TestUtils.containsAll(err.toString(),
          "Skipping record with no RG id",
          "Skipped 6 records with no RG id"
          );

      final ReadGroupStatsCalculator calculator2 = new ReadGroupStatsCalculator();
      calculator2.merge(calculator);
      calculator2.calculate();
      final ReadGroupStats rgs1 = calculator.getStats("boo_pe");
      final ReadGroupStats rgs2 = calculator2.getStats("boo_pe");

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertEquals(rgs1.fragmentStdDev(), rgs2.fragmentStdDev());

      calculator2.merge(calculator);
      calculator2.calculate();

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertTrue(rgs1.fragmentStdDev() > rgs2.fragmentStdDev());
      assertEquals(41.79664, rgs2.fragmentStdDev(), 0.00001);

      calculator2.merge(calculator);
      calculator2.calculate();

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertTrue(rgs1.fragmentStdDev() > rgs2.fragmentStdDev());
      assertEquals(41.59901, rgs2.fragmentStdDev(), 0.00001);
    }
  }

  public void testBlend() throws IOException {
    final ReadGroupStatsCalculator.Merger merger = new ReadGroupStatsCalculator.Merger();
    try (final TestDirectory tempDir = new TestDirectory("readgroup")) {
      final File samFile = FileUtils.stringToFile(RG_SAM, new File(tempDir, "rg.sam"));
      final List<File> files = new ArrayList<>();
      files.add(samFile);
      final ReadGroupStatsCalculator calculator = merger.createReadGroupStatsCalculator();
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      calculator.calculate(files, out.outputStream());
      //System.err.println(out.toString());
      //System.err.println(err.toString());
      TestUtils.containsAll(out.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "boo_pe\t"
          );
      TestUtils.containsAll(err.toString(),
          "Skipping record with no RG id",
          "Skipped 6 records with no RG id"
          );

      ReadGroupStatsCalculator calcBlend = merger.blend();
      final ReadGroupStats rgs1 = calculator.getStats("boo_pe");
      ReadGroupStats rgs2 = calcBlend.getStats("boo_pe");

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertEquals(rgs1.fragmentStdDev(), rgs2.fragmentStdDev());

      merger.createReadGroupStatsCalculator().calculate(files, out.outputStream());

      calcBlend = merger.blend();
      rgs2 = calcBlend.getStats("boo_pe");

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertTrue(rgs1.fragmentStdDev() > rgs2.fragmentStdDev());
      assertEquals(41.79664, rgs2.fragmentStdDev(), 0.00001);

      merger.createReadGroupStatsCalculator().calculate(files, out.outputStream());

      calcBlend = merger.blend();
      rgs2 = calcBlend.getStats("boo_pe");

      assertEquals(rgs1.fragmentMean(), rgs2.fragmentMean());
      assertTrue(rgs1.fragmentStdDev() > rgs2.fragmentStdDev());
      assertEquals(41.59901, rgs2.fragmentStdDev(), 0.00001);
    }
  }
}
