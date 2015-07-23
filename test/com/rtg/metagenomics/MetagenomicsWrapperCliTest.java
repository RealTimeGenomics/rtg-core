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
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.LoggedCli;
import com.rtg.metagenomics.MetagenomicsWrapperCli.Platform;
import com.rtg.ngs.MapFlags;
import com.rtg.protein.MapXCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.LogStream;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

/**
 */
public class MetagenomicsWrapperCliTest extends AbstractCliTest {

  public void testPlatformEnum() {
    TestUtils.testEnum(Platform.class, "[ILLUMINA, IONTORRENT]");
  }

  public void testFlags() throws IOException {
    try (TestDirectory tmp = new TestDirectory()) {
      final File protein = new File(tmp, "protein");
      final File dna = new File(tmp, "dna");
      final File readSdf = new File(tmp, "dna");
      final File iFile = new File(tmp, "file");
      assertTrue(iFile.createNewFile());
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGA" + LS, dna, new SdfId());

      String err = checkHandleFlagsErr("--protein", protein.getPath(), "--species", dna.getPath(), "--input", iFile.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Filter reference required (use --filter)");
      err = checkHandleFlagsErr("--protein", protein.getPath(), "--filter", dna.getPath(), "--input", iFile.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Species reference required (use --species)");
      err = checkHandleFlagsErr("--species", dna.getPath(), "--filter", dna.getPath(), "--input", iFile.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Protein reference required (use --protein)");
      err = checkHandleFlagsErr("--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Some read data is required");
      err = checkHandleFlagsErr("--input-left", "LEFT", "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "You must provide both --input-left and --input-right");
      err = checkHandleFlagsErr("--input-left", iFile.getPath(), "--input-right", "RIGHT", "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "--input-right should be a file");
      err = checkHandleFlagsErr("--input-left", "LEFT", "--input-right", iFile.getPath(), "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "--input-left should be a file");
      err = checkHandleFlagsErr("--input", iFile.getPath(), "--input-right", iFile.getPath(), "--input-left", iFile.getPath(), "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Too many read datasets supplied.");
      err = checkHandleFlagsErr("--input", "FILE", "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "--input file doesn't exist: FILE");
      err = checkHandleFlagsErr("--input", readSdf.getPath(), "--protein", "PROT", "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "The specified SDF, \"PROT\", does not exist.");
      err = checkHandleFlagsErr("--input", readSdf.getPath(), "--protein", protein.getPath(), "--species", "SPEC", "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "The specified SDF, \"SPEC\", does not exist.");
      err = checkHandleFlagsErr("--input", readSdf.getPath(), "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", "FILT", "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "The specified SDF, \"FILT\", does not exist.");
      err = checkHandleFlagsErr("--input", readSdf.getPath(), "--protein", tmp.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Invalid protein database.");
      err = checkHandleFlagsErr("--input", readSdf.getPath(), "--protein", dna.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
      TestUtils.containsAll(err, "Protein database should be an SDF containing formatted protein data");
      checkHandleFlagsOut("--input", readSdf.getPath(), "--protein", protein.getPath(), "--species", dna.getPath(), "--filter", dna.getPath(), "--output", new File(tmp, "output").getPath());
    }
  }

  public void testSimple() throws IOException {
    try (TestDirectory tmp = new TestDirectory()) {
      final File protein = new File(tmp, "protein");
      final File dna = new File(tmp, "dna");
      final File filter = new File(tmp, "filter");
      final File readSdf = new File(tmp, "reads");
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, dna, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, readSdf, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      final File output = new File(tmp, "output");
      final String[] flags = {"--input", readSdf.getPath(), "--filter", filter.getPath(), "--species", dna.getPath(), "--protein", protein.getPath(), "--output", output.getPath()};
      final String res = checkMainInitOk(flags);
      TestUtils.containsAll(res
          , "## rtg mapf "
          , "## rtg map "
          , "## rtg species "
          , "## rtg mapx "
          );
    }
  }

  public void testDefaults() throws IOException {
    final String oldRefDir = System.getProperty("references.dir");
    try (TestDirectory tmp = new TestDirectory()) {
      System.setProperty("references.dir", tmp.getPath());
        final CFlags flags = new CFlags();
        MetagenomicsWrapperCli.initFlags(flags, true, true);
        assertEquals(new File(tmp, "filter").getPath(), flags.getFlag("filter").getParameterDefault().toString());
        assertEquals(new File(tmp, "species").getPath(), flags.getFlag("species").getParameterDefault().toString());
        assertEquals(new File(tmp, "protein").getPath(), flags.getFlag("protein").getParameterDefault().toString());
    } finally {
      if (oldRefDir != null) {
        System.setProperty("references.dir", oldRefDir);
      } else {
        System.clearProperty("references.dir");
      }
    }
  }

  public void testGeneratedFlags() throws IOException {
    try (TestDirectory tmp = new TestDirectory()) {
      final File protein = new File(tmp, "protein");
      final File dna = new File(tmp, "dna");
      final File filter = new File(tmp, "filter");
      final File readSdf = new File(tmp, "reads");
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGTACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, dna, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGACAGTTTAGGAAAAAAAAAAAATTTTTTTTTTTTTTTCCCCCCCCCGGGGGGGGGGGATGCACGT" + LS, readSdf, new SdfId());
      ReaderTestUtils.getReaderDNA(">a" + LS + "AAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCCAAAATAAAGGAAAGGTTTCC" + LS, filter, new SdfId());
      final File output = new File(tmp, "output");
      final String[] flags = {"--input", readSdf.getPath(), "--filter", filter.getPath(), "--species", dna.getPath(), "--protein", protein.getPath(), "--output", output.getPath()};


      final List<List<String>> expected = new ArrayList<>();
      final File mapf = new File(output, "mapf");
      final File map = new File(output, "map");
      final File species = new File(output, "species");
      final File mapx = new File(output, "mapx1");
      expected.add(Arrays.asList("mapf", "--output", mapf.getPath(), "--template", filter.getPath(), "--sam-rg", "@RG\\tPL:ILLUMINA\\tSM:sample\\tID:id", "--input", readSdf.getPath()));
      expected.add(Arrays.asList("map", "--output", map.getPath(), "--template", dna.getPath(), "--" + MapFlags.MAX_ALIGNMENT_MISMATCHES, "10%", "--max-top-results", "100", "--sam-rg", "@RG\\tPL:ILLUMINA\\tSM:sample\\tID:id", "--input", new File(mapf, "unmapped.sdf").getPath()));
      expected.add(Arrays.asList("species", "--output", species.getPath(), "--genomes", dna.getPath(), new File(map, "bazinga.bam").getPath(), new File(map, "random.bam").getPath()));
      expected.add(Arrays.asList("mapx", "--output", mapx.getPath(), "--template", protein.getPath(), "--" + MapXCli.MAX_ALIGNMENT_SCORE, "10%", "--max-top-results", "10", "--input", new File(mapf, "unmapped.sdf").getPath(), "--format", "fastq"));

      final List<List<String>> commands = generateCommands(flags, output);
      final List<String> modules = new ArrayList<>();
      for (List<String> command : commands) {
        modules.add(command.get(0));
      }
      assertEquals(modules.toString(), expected.size(), commands.size());
      for (int i = 0; i < expected.size(); i++) {
        assertEquals(expected.get(i), commands.get(i));

      }
    }
  }

  private static class MockMetagenomics extends MetagenomicsWrapperCli {

    private final File mOutputDir;
    private final List<List<String>> mCommands;

    public MockMetagenomics(File outputDir, List<List<String>> commands) {
      mOutputDir = outputDir;
      mCommands = commands;
    }

    private class MockMetagenomicsWrapperTask extends MetagenomicsWrapperTask {
      protected MockMetagenomicsWrapperTask(MetaPipelineParams params, OutputStream reportStream, UsageMetric usageMetric, LogStream logStream, PrintStream err) {
        super(params, reportStream, usageMetric, logStream, err);
      }

      @Override
      public int runCommand(OutputStream out, LoggedCli module, List<String> args) throws IOException {
        // Create mapf output so that a file can be found for mapping/protein
        // Created here to prevent the wrapper from bailing due to non empty directories
        final File mapfDir = new File(mOutputDir, "map");
        final File sdfFile1 = new File(mapfDir, "random.bam");
        final File sdfFile2 = new File(mapfDir, "bazinga.bam");
        if (!mapfDir.exists()) {
          boolean success = true;
          success &= mapfDir.mkdir();
          success &= sdfFile1.createNewFile();
          success &= sdfFile2.createNewFile();
          if (!success) {
            throw new IOException("Can't create mock files");
          }
        }
        final List<String> flags = new ArrayList<>();
        flags.add(module.moduleName());
        flags.addAll(args);

        mCommands.add(flags);
        return 0;
      }
    }

    @Override
    protected IORunnable task(MetaPipelineParams params, OutputStream out) {
      assertFalse(params.outputParams().isCompressed());
      assertFalse(params.outputParams().progress());
      mTask = new MockMetagenomicsWrapperTask(params, out, mUsageMetric, mLogStream, mErr);
      return mTask;
    }
  }

  private List<List<String>> generateCommands(String[] flags, File output) {
    final List<List<String>> commands = new ArrayList<>();
    final MetagenomicsWrapperCli cli = new MockMetagenomics(output, commands);
    final MemoryPrintStream mps = new MemoryPrintStream();
    cli.mainInit(flags, mps.outputStream(), mps.printStream());
    return commands;
  }

  @Override
  protected AbstractCli getCli() {
    return new MetagenomicsWrapperCli();
  }
}
