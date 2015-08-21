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

package com.rtg.metagenomics;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.metagenomics.MetagenomicsWrapperCli.Platform;
import com.rtg.ngs.MapFlags;
import com.rtg.protein.MapXCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.usage.UsageMetric;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogSimple;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class MetagenomicsWrapperTaskTest extends TestCase {

  private CliDiagnosticListener mListener;
  private MemoryPrintStream mOut;
  private MemoryPrintStream mErr;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
    CommandLine.clearCommandArgs();
    mOut = new MemoryPrintStream();
    mErr = new MemoryPrintStream();
    mListener = new CliDiagnosticListener(mErr.printStream(), mOut.printStream());
    Diagnostic.addListener(mListener);
  }

  @Override
  protected void tearDown() {
    Diagnostic.removeListener(mListener);
    mOut = null;
    mErr = null;
    mListener = null;
  }

  public void testErrorRate() {
    assertEquals("10%", MetagenomicsWrapperTask.errorRate(Platform.ILLUMINA));
    assertEquals("15%", MetagenomicsWrapperTask.errorRate(Platform.IONTORRENT));
  }

  public void testBasic1() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File filter = new File(tmp, "filter");
      final File species = new File(tmp, "species");
      final File protein = new File(tmp, "protein");
      final File readSdf = new File(tmp, "reads");
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, species, new SdfId());
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, readSdf, new SdfId());
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .filterSdf(filter)
          .speciesSdf(species)
          .proteinSdf(protein)
          .inputFile(readSdf)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(mErr.toString(), 0, task.returnCode());
      task.run();
      assertEquals(mErr.toString(), 0, task.returnCode());
      final File mapfOut = new File(output, "mapf");
      final File mapxOut = new File(output, "mapx1");
      final File mapOut = new File(output, "map");
      final File speciesOut = new File(output, "species");
      final File reportOut = new File(output, "report");
      assertTrue(mapfOut.isDirectory());
      assertTrue(mapxOut.isDirectory());
      assertTrue(mapOut.isDirectory());
      assertTrue(speciesOut.isDirectory());
      assertTrue(reportOut.isDirectory());
      TestUtils.containsAll(mOut.toString()
          , "## rtg mapf  --output " + mapfOut.getPath() + " --template " + filter.getPath() + " --sam-rg @RG\\tPL:ILLUMINA\\tSM:sample\\tID:id --input " + readSdf.getPath()
          , "## rtg map  --output " + mapOut.getPath() + " --template " + species.getPath() + " --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " 10% --max-top-results 100 --sam-rg @RG\\tPL:ILLUMINA\\tSM:sample\\tID:id --input " + new File(mapfOut, "unmapped.sdf").getPath()
          , "## rtg species  --output " + speciesOut.getPath() + " --genomes " + species.getPath() + " " + new File(mapOut, "alignments.bam").getPath()
          , "## rtg mapx  --output " + mapxOut.getPath() + " --template " + protein.getPath() + " --" + MapXCli.MAX_ALIGNMENT_SCORE + " 10% --max-top-results 10 --input " + new File(mapfOut, "unmapped.sdf").getPath()
          );
    }
  }

  public void testBasic2() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File filter = new File(tmp, "filter");
      final File species = new File(tmp, "species");
      final File protein = new File(tmp, "protein");
      final File readsLeft = new File(tmp, "reads-left.fastq");
      final File readsRight = new File(tmp, "reads-right.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, species, new SdfId());
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
                           + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, readsLeft);
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
                           + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, readsRight);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .filterSdf(filter)
          .speciesSdf(species)
          .proteinSdf(protein)
          .inputLeft(readsLeft)
          .inputRight(readsRight)
          .inputPlatform(Platform.IONTORRENT)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(mErr.toString(), 0, task.returnCode());
      final File mapfOut = new File(output, "mapf");
      final File mapxOut1 = new File(output, "mapx1");
      final File mapxOut2 = new File(output, "mapx2");
      final File mapOut = new File(output, "map");
      final File speciesOut = new File(output, "species");
      final File reportOut = new File(output, "report");
      assertTrue(mapfOut.isDirectory());
      assertTrue(mapxOut1.isDirectory());
      assertTrue(mapxOut2.isDirectory());
      assertTrue(mapOut.isDirectory());
      assertTrue(speciesOut.isDirectory());
      assertTrue(reportOut.isDirectory());
      final File unmappedSdf = new File(mapfOut, "unmapped.sdf");
      TestUtils.containsAll(mOut.toString()
          , "## rtg mapf  --output " + mapfOut.getPath() + " --template " + filter.getPath() + " --sam-rg @RG\\tPL:IONTORRENT\\tSM:sample\\tID:id --format fastq --quality-format sanger --left " + readsLeft.getPath() + " --right " + readsRight.getPath()
          , "## rtg map  --output " + mapOut.getPath() + " --template " + species.getPath() + " --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " 15% --max-top-results 100 --sam-rg @RG\\tPL:IONTORRENT\\tSM:sample\\tID:id --input " + unmappedSdf.getPath()
          , "## rtg species  --output " + speciesOut.getPath() + " --genomes " + species.getPath() + " " + new File(mapOut, "alignments.bam").getPath()
          , "## rtg mapx  --output " + mapxOut1.getPath() + " --template " + protein.getPath() + " --" + MapXCli.MAX_ALIGNMENT_SCORE + " 15% --max-top-results 10 --input " + new File(unmappedSdf, "left").getPath()
          , "## rtg mapx  --output " + mapxOut2.getPath() + " --template " + protein.getPath() + " --" + MapXCli.MAX_ALIGNMENT_SCORE + " 15% --max-top-results 10 --input " + new File(unmappedSdf, "right").getPath()
          );
    }
  }

  public void testNoFilter() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File species = new File(tmp, "species");
      final File protein = new File(tmp, "protein");
      final File reads = new File(tmp, "reads");
      assertTrue(reads.mkdir());
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, species, new SdfId());
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.createPairedReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS,  ">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, reads, new SdfId());
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .speciesSdf(species)
          .proteinSdf(protein)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(mErr.toString(), 0, task.returnCode());
      final File mapOut = new File(output, "map");
      final File speciesOut = new File(output, "species");
      final File mapxOut1 = new File(output, "mapx1");
      final File mapxOut2 = new File(output, "mapx2");
      final File reportOut = new File(output, "report");
      assertTrue(mapOut.isDirectory());
      assertTrue(speciesOut.isDirectory());
      assertTrue(mapxOut1.isDirectory());
      assertTrue(mapxOut2.isDirectory());
      assertTrue(reportOut.isDirectory());
      TestUtils.containsAll(mOut.toString()
          , "## rtg map  --output " + mapOut.getPath() + " --template " + species.getPath() + " --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " 10% --max-top-results 100 --sam-rg @RG\\tPL:ILLUMINA\\tSM:sample\\tID:id --input " + reads.getPath()
          , "## rtg species  --output " + speciesOut.getPath() + " --genomes " + species.getPath() + " " + new File(mapOut, "alignments.bam").getPath()
          , "## rtg mapx  --output " + mapxOut1.getPath() + " --template " + protein.getPath() + " --" + MapXCli.MAX_ALIGNMENT_SCORE + " 10% --max-top-results 10 --input " + new File(reads, "left").getPath()
          , "## rtg mapx  --output " + mapxOut2.getPath() + " --template " + protein.getPath() + " --" + MapXCli.MAX_ALIGNMENT_SCORE + " 10% --max-top-results 10 --input " + new File(reads, "right").getPath()
          );
    }
  }

  public void testFilterInputFlags1() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File filter = new File(tmp, "filter");
      final File readsLeft = new File(tmp, "reads-left.fastq");
      final File readsRight = new File(tmp, "reads-right.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
                           + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, readsLeft);
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
                           + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, readsRight);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .filterSdf(filter)
          .inputLeft(readsLeft)
          .inputRight(readsRight)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(0, task.returnCode());
      final File mapfOut = new File(output, "mapf");
      final File reportOut = new File(output, "report");
      assertTrue(mapfOut.isDirectory());
      assertTrue(reportOut.isDirectory());
      TestUtils.containsAll(mOut.toString()
          , "## rtg mapf  --output " + mapfOut.getPath() + " --template " + filter.getPath() + " --sam-rg @RG\\tPL:ILLUMINA\\tSM:sample\\tID:id --format fastq --quality-format sanger --left " + readsLeft.getPath() + " --right " + readsRight.getPath()
          );
    }
  }

  public void testFilterInputFlags2() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File filter = new File(tmp, "filter");
      final File reads = new File(tmp, "reads.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
                           + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, reads);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .filterSdf(filter)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(0, task.returnCode());
      final File mapfOut = new File(output, "mapf");
      final File reportOut = new File(output, "report");
      assertTrue(mapfOut.isDirectory());
      assertTrue(reportOut.isDirectory());
      TestUtils.containsAll(mOut.toString()
          , "## rtg mapf  --output " + mapfOut.getPath() + " --template " + filter.getPath() + " --sam-rg @RG\\tPL:ILLUMINA\\tSM:sample\\tID:id --format fastq --quality-format sanger --input " + reads
          );
    }
  }

  public void testFilterFail() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File filter = new File(tmp, "filter");
      final File reads = new File(tmp, "reads.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      FileUtils.stringToFile(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, reads);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .filterSdf(filter)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(1, task.returnCode());
      final File reportOut = new File(output, "report");
      assertFalse(reportOut.isDirectory());
      TestUtils.containsAll(mErr.toString(), "At least one input file looks like FASTA");
    }
  }

  public void testMapFail() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File species = new File(tmp, "species");
      final File reads = new File(tmp, "reads.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, species, new SdfId());
      FileUtils.stringToFile(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, reads);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .speciesSdf(species)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(1, task.returnCode());
      final File reportOut = new File(output, "report");
      assertFalse(reportOut.isDirectory());
      TestUtils.containsAll(mErr.toString(), "At least one input file looks like FASTA");
    }
  }

  public void testSpeciesFail() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File species = new File(tmp, "species");
      final File reads = new File(tmp, "reads.fastq");
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, species, new SdfId());
      final File taxonFile = new File(species, "taxonomy.tsv");
      assertTrue(taxonFile.createNewFile());
      FileUtils.stringToFile("@a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS
          + "+a" + LS + "2222222222222222222222222222222222222222222222222222222222222222222222222" + LS, reads);
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .speciesSdf(species)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(1, task.returnCode());
      final File mapOut = new File(output, "map");
      final File reportOut = new File(output, "report");
      assertTrue(mapOut.isDirectory());
      assertFalse(reportOut.isDirectory());
      TestUtils.containsAll(mErr.toString(), "Reference SDF does not contain both taxonomy and sequences lookup");
    }
  }

  public void testMapxFail() throws IOException {
    try (final TestDirectory tmp = new TestDirectory()) {
      final File protein = new File(tmp, "protein");
      final File reads = new File(tmp, "reads");
      assertTrue(reads.mkdir());
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.createPairedReaderDNA(">a" + LS + "ACGTTTAGACA" + LS,  ">a" + LS + "ACGTTTAGA" + LS, reads, new SdfId());
      final File output = new File(tmp, "output");
      final MetaPipelineParams params = MetaPipelineParams.builder()
          .proteinSdf(protein)
          .inputFile(reads)
          .outputParams(new OutputParams(output, false, false))
          .create();
      final MemoryPrintStream logStream = new MemoryPrintStream();
      final LogSimple log = new LogSimple(logStream.printStream());
      final MetagenomicsWrapperTask task = new MetagenomicsWrapperTask(params, mOut.outputStream(), new UsageMetric(), log, mErr.printStream());
      assertEquals(0, task.returnCode());
      task.run();
      assertEquals(1, task.returnCode());
      final File reportOut = new File(output, "report");
      assertFalse(reportOut.isDirectory());
      TestUtils.containsAll(mErr.toString(), "The word length \"7\" should be less than the read length \"2\"");
    }
  }
}
