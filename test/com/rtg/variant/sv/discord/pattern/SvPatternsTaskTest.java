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
import java.util.Collections;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.header.VcfHeader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

public class SvPatternsTaskTest extends AbstractNanoTest {

  public void testEnd2End() throws IOException, UnindexableDataException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final SAMFileHeader samHeader = new SAMFileHeader();
      samHeader.addSequence(new SAMSequenceRecord("s1", 1000));
      samHeader.addSequence(new SAMSequenceRecord("s2", 1000));
      final VcfHeader vcfHeader = new VcfHeader();
      vcfHeader.addContigFields(samHeader);

      final File vcfFile = new File(tmpDir, "input.vcf.gz");
      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", vcfFile);
      final TabixIndexer tabixIndexer = new TabixIndexer(vcfFile, new File(vcfFile.getPath() + TabixIndexer.TABIX_EXTENSION));
      tabixIndexer.saveVcfIndex();
      final File output = new File(tmpDir, "output");
      assertTrue(output.mkdir());
      CommandLine.setCommandArgs("foo", "bar", "baz");
      final BreakpointPatternParams params = BreakpointPatternParams.builder().directory(output).files(Collections.singletonList(vcfFile)).create();
      final MemoryPrintStream mps = new MemoryPrintStream();
      final SvPatternsTask task = new SvPatternsTask(params, mps.outputStream());
      task.exec();
      assertTrue(output.isDirectory());
      final File bedFile = new File(output, "sv_patterns.bed.gz");
      assertTrue(bedFile.isFile());
      assertTrue(new File(output, "sv_patterns.bed.gz.tbi").isFile());
      final String results = FileHelper.gzFileToString(bedFile);
      TestUtils.containsAll(results, "SV Patterns output 0.1");
      mNano.check("sv_patterns.bed", results.replaceAll("#Version.*(\r)?\n", "#Version[...]\n"));
    }
  }

  public void testLoad() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File vcfFile = new File(tmpDir, "input.vcf.gz");
      // This file contains a few records that ought to be filtered as well as the ones from {@code EXPECTED_BREAKPOINTS}
      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", vcfFile);
      final BreakpointPatternParams params = BreakpointPatternParams.builder().files(Collections.singletonList(vcfFile)).create();
      final MemoryPrintStream mps = new MemoryPrintStream();
      final SvPatternsTask task = new SvPatternsTask(params, mps.outputStream());
      final BreakpointStore store = task.loadBreakpoints(null, Collections.singletonList(vcfFile));
      final StringBuilder sb = new StringBuilder();
      for (VcfBreakpoint b : store) {
        sb.append(b.toString()).append(LS);
      }
      mNano.check("sv_load.txt", sb.toString());

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  String tab(String subject) {
    return subject.replaceAll(" ", "\t");
  }

  public void testAnalyzeDelete() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, false, 0));
    final String results = processStore(store);
    assertEquals(tab("a 11 5000 deletion 0" + LS), results);
  }

  public void testOverlapDelete() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, false, 0));
    store.add(new VcfBreakpoint("a", 400, "a", 5000, true, false, 2));
    final String results = processStore(store);
    assertEquals(tab("a 401 5000 deletion 2" + LS), results);

  }

  public void testAnalyzeInversion() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 1));
    store.add(new VcfBreakpoint("a", 14, "a", 5002, false, false, 1));
    final String results = processStore(store);
    assertEquals(tab("a 10 5000 inversion 2" + LS), results);
  }
  public void testSlightlyOffInversion() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "a", 5002, false, false, 1));
    final String results = processStore(store);
    assertEquals(tab("a 10 5000 inversion 1" + LS), results);
  }
  public void testNotAnInversion() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "a", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "a", 5002, false, true, 0));
    final String results = processStore(store);
    assertEquals("", results);
  }


  public void testAnalyzeCopy() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 3));
    store.add(new VcfBreakpoint("a", 14, "b", 4700, false, false, 0));
    final String results = processStore(store);
    assertEquals(tab("a 10 14 inserted_copy:b:5000-4700 3" + LS), results);
  }

  public void testSlightlyOffCopy() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 8, "b", 4700, false, false, 0));
    final String results = processStore(store);
    assertEquals(tab("a 8 10 inserted_copy:b:5000-4700 0" + LS), results);
  }
  public void testNotACopy() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 5000, true, true, 0));
    store.add(new VcfBreakpoint("a", 14, "b", 4700, false, true, 0));
    final String results = processStore(store);
    assertEquals("", results);
  }

  public void testForwardCopy() throws IOException {
    final BreakpointStore store = new BreakpointStore();
    store.add(new VcfBreakpoint("a", 10, "b", 4700, true, false, 1));
    store.add(new VcfBreakpoint("a", 14, "b", 5000, false, true, 5));
    final String results = processStore(store);
    assertEquals(tab("a 10 14 inserted_copy:b:4700-5000 6" + LS), results);
  }

  public String processStore(BreakpointStore store) throws IOException {
    final BreakpointPatternParams params = BreakpointPatternParams.builder().create();
    final MemoryPrintStream stdOut = new MemoryPrintStream();
    final MemoryPrintStream bedFile = new MemoryPrintStream();
    final SvPatternsTask task = new SvPatternsTask(params, stdOut.outputStream());
    task.setOutput(bedFile.outputStream());
    for (VcfBreakpoint breakpoint : store) {
      task.analyseBreakpoint(breakpoint, store);
    }
    task.close();
    return bedFile.toString();
  }

}
