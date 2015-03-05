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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.header.VcfHeader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 *         Date: 19/03/12
 *         Time: 10:10 AM
 */
public class SvPatternsTaskTest extends TestCase {

  public void testEnd2End() throws IOException, UnindexableDataException {
    CommandLine.setCommandArgs("foo", "bar", "baz");
    try {
      final SAMFileHeader samHeader = new SAMFileHeader();
      samHeader.addSequence(new SAMSequenceRecord("s1", 1000));
      samHeader.addSequence(new SAMSequenceRecord("s2", 1000));
      final VcfHeader vcfHeader = new VcfHeader();
      vcfHeader.addContigFields(samHeader);

      final File tmpDir = FileHelper.createTempDirectory();
      try {
        final File vcfFile = new File(tmpDir, "input.vcf.gz");
        FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", vcfFile);
        TabixIndexer tabixIndexer = new TabixIndexer(vcfFile, new File(vcfFile.getPath() + TabixIndexer.TABIX_EXTENSION));
        tabixIndexer.saveVcfIndex();
        final File output = new File(tmpDir, "output");
        assertTrue(output.mkdir());
        final BreakpointPatternParams params = BreakpointPatternParams.builder().directory(output).files(Arrays.asList(vcfFile)).create();
        MemoryPrintStream mps = new MemoryPrintStream();
        SvPatternsTask task = new SvPatternsTask(params, mps.outputStream());
        task.exec();
        assertTrue(output.isDirectory());
        final File bedFile = new File(output, "sv_patterns.bed.gz");
        assertTrue(bedFile.isFile());
        assertTrue(new File(output, "sv_patterns.bed.gz.tbi").isFile());
        String results = FileHelper.gzFileToString(bedFile);
        TestUtils.containsAll(results,
            "#CL\tfoo bar baz"
            , "#Version"
            , "SV Patterns output 0.1"
            , "#sequence\tstart\tend\tdescription\tcount" + LS
        );


      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    } finally {
      CommandLine.clearCommandArgs();
    }
  }

  static final String[] EXPECTED_BREAKPOINTS = {
      "VcfBreakpoint: simulatedSequence1 50000 simulatedSequence1 53016 true true"
      , "VcfBreakpoint: simulatedSequence1 50003 simulatedSequence1 53013 false false"
      , "VcfBreakpoint: simulatedSequence1 53016 simulatedSequence1 50000 true true"
      , "VcfBreakpoint: simulatedSequence1 53013 simulatedSequence2 50003 false false"
      , "VcfBreakpoint: simulatedSequence2 53016 simulatedSequence1 50000 true true"
      , "VcfBreakpoint: simulatedSequence2 53015 simulatedSequence2 50000 true true"
  };
  public void testLoad() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File vcfFile = new File(tmpDir, "input.vcf.gz");
      // This file contains a few records that ought to be filtered as well as the ones from {@code EXPECTED_BREAKPOINTS}
      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", vcfFile);
      final BreakpointPatternParams params = BreakpointPatternParams.builder().files(Arrays.asList(vcfFile)).create();
      MemoryPrintStream mps = new MemoryPrintStream();
      SvPatternsTask task = new SvPatternsTask(params, mps.outputStream());
      BreakpointStore store = task.loadBreakpoints(null, Arrays.asList(vcfFile));
      int i = 0;
      for (VcfBreakpoint b : store) {
        assertEquals(EXPECTED_BREAKPOINTS[i], b.toString());
        i++;
      }
      assertEquals(6, i);

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  String tab(String subject) {
    return subject.replaceAll(" ", "\t");
  }

  public void testAnalyzeDelete() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, false, 0));
    String results = processStore(store);
    assertEquals(tab("a 10 5000 deletion 0" + LS), results);
  }

  public void testOverlapDelete() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, false, 0));
    store.add(new VcfBreakpoint("a", 400, "a", 5000, true, false, 2));
    String results = processStore(store);
    assertEquals(tab("a 400 5000 deletion 2" + LS), results);

  }

  public void testAnalyzeInversion() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 1));
    store.add(new VcfBreakpoint("a", 14, "a", 5002, false, false, 1));
    String results = processStore(store);
    assertEquals(tab("a 10 5000 inversion 2" + LS), results);
  }
  public void testSlightlyOffInversion() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "a", 5002, false, false, 1));
    String results = processStore(store);
    assertEquals(tab("a 10 5000 inversion 1" + LS), results);
  }
  public void testNotAnInversion() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "a", 5002, false, true, 0));
    String results = processStore(store);
    assertEquals("", results);
  }


  public void testAnalyzeCopy() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 3));
    store.add(new VcfBreakpoint("a", 14, "b", 4700, false, false, 0));
    String results = processStore(store);
    assertEquals(tab("a 10 14 inserted_copy:b:5000-4700 3" + LS), results);
  }

  public void testSlightlyOffCopy() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "b", 4700, false, false, 0));
    String results = processStore(store);
    assertEquals(tab("a 8 10 inserted_copy:b:5000-4700 0" + LS), results);
  }
  public void testNotACopy() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 14, "b", 4700, false, true, 0));
    String results = processStore(store);
    assertEquals("", results);
  }

  public void testForwardCopy() throws IOException {
    BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 4700, true, false, 1));
    store.add(new VcfBreakpoint("a", 14, "b", 5000, false, true, 5));
    String results = processStore(store);
    assertEquals(tab("a 10 14 inserted_copy:b:4700-5000 6" + LS), results);
  }

  public String processStore(BreakpointStore store) throws IOException {
    final BreakpointPatternParams params = BreakpointPatternParams.builder().create();
    MemoryPrintStream stdOut = new MemoryPrintStream();
    MemoryPrintStream bedFile = new MemoryPrintStream();
    SvPatternsTask task = new SvPatternsTask(params, stdOut.outputStream());
    task.setOutput(bedFile.outputStream());
    task.processStore(store);
    task.close();
    return bedFile.toString();
  }

}
