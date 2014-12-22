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
package com.rtg.simulation.snpsim;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.BuildTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;

/**
 * Test the corresponding class
 */
public class GenomeMutatorCliTest extends AbstractCliTest {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  @Override
  protected AbstractCli getCli() {
    return getGMCli();
  }

  private GenomeMutatorCli getGMCli() {
    return new GenomeMutatorCli() {
      @Override
      GenomeMutator getMutator(Integer seed, boolean verbose, boolean simpleMnps, int minDistance, VariantParams params, String sample) {
        return new GenomeMutator(verbose, simpleMnps, 1, params, sample) {
          @Override
          int mutatePriors(SequencesReader input, Sex sex, File outputDirectory, File twinDirectory, OutputStream mappingOutput) throws IOException {
            if (twinDirectory != null) {
              return -4;
            }
            mFormatter.writeHeader(mappingOutput, mVariantParams, null);
            return -1;
          }
          @Override
          int mutateCount(final SequencesReader input, Sex sex, final File outputDirectory, final File twinDirectory, final OutputStream mappingOutput, final int count) {
            return -2;
          }
          @Override
          int mutateRate(final SequencesReader input, Sex sex, final File outputDirectory, final File twinDirectory, final OutputStream mappingOutput, final double rate) {
            assertEquals(0.5, rate);
            return -3;
          }
          @Override
          void setRates(final Double snpRate, final Double mnpRate, final Double indelRate) {
            assertEquals(0.0, snpRate);
            assertEquals(0.0, mnpRate);
            assertEquals(0.0, indelRate);
            super.setRates(snpRate, mnpRate, indelRate);
          }
          @Override
          void summarise(File outputDir, File twinOutputDir, OutputStream out) {
          }
        };
      }
    };
  }

  private static final OutputStream NULL_STREAM = TestUtils.getNullOutputStream();
  //private static final PrintStream NULL_PRINTSTREAM = SimpleTestUtils.getNullPrintStream();

  public void testName() {
    final GenomeMutatorCli sim = getGMCli();
    assertEquals("rtg snpsim", sim.applicationName() + " " + sim.moduleName());
  }

  /**
   * Test method for {@link GenomeMutatorCli}.
   */
  public final void testGetCFlags() {
    final CFlags flags = new CFlags();
    final GenomeMutatorCli sim = getGMCli();
    sim.initFlags(flags);

    TestCFlags.check(flags, //"rtg snpsim [OPTION]... -i SDF -o SDF -s FILE",
        "input=", "input genome",
        "-o", "output=", "output SDF",
        "snp-file", "output file with SNP information",
        "print help on command-line flag usage",
        //"selects a properties file specifying the priors. Either a",
        //"file name or one of [human] (Default is human)",
        "seed=", "seed for the random number generator",
        "-Z", "no-gzip",
        "-O", "--diploid", "secondary genome");

        TestCFlags.checkExtendedUsage(flags,
        "fraction of nucleotide positions to be altered (default",
        "number of mutations to be inserted",
        "fraction of mutations that are SNP",
        "fraction of mutations that are MNP",
        "output into SNP file used for debugging");
  }

  public void testSomeErrors() throws IOException {
    final File tmp = FileUtils.createTempDir("GenomeMutator", "tmp");
    try {
      final File mutant = new File(tmp, "mutant");
      final File mutations = new File(tmp, "mutations");
      final File in = BuildTestUtils.prereadDNA(mDir, ">a\nacgtacgatcagcatctgac\n");

      final GenomeMutatorCli cli = getGMCli();
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      final PrintStream errorPrint = new PrintStream(baos);
      assertEquals(1, cli.mainInit(new String[] {"-o", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath(), "-Z", "-p", "blob"}, NULL_STREAM, errorPrint));
      errorPrint.flush();
      assertTrue(baos.toString().contains("Invalid prior option \"blob\""));
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testFlags() throws IOException {
    final MemoryPrintStream errorPrint = new MemoryPrintStream();
    Diagnostic.setLogStream(errorPrint.printStream());

    final File tmp = FileUtils.createTempDir("GenomeMutator", "tmp");
    try {
      final File mutant = new File(tmp, "mutant");
      final File mutations = new File(tmp, "mutations");
      final File in = BuildTestUtils.prereadDNA(mDir, ">a\nacgtacgatcagcatctgac\n");

      final GenomeMutatorCli cli = getGMCli();

      final ByteArrayOutputStream baos = new ByteArrayOutputStream();

      //test use priors
      int outCode = cli.mainInit(new String[] {"-o", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath(), "-Z"}, NULL_STREAM, errorPrint.printStream());
      assertEquals(errorPrint.toString(), -1, outCode);
      errorPrint.reset();
      assertTrue(mutations.exists());
      assertTrue(!mutant.exists() || FileHelper.deleteAll(mutant));
      //test gzip output
      outCode = cli.mainInit(new String[] {"-o", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath()}, errorPrint.printStream(), errorPrint.printStream());
      assertEquals(errorPrint.toString(), -1, outCode);
      assertTrue(new File(mutations.getAbsolutePath() + FileUtils.GZ_SUFFIX).exists());
      assertTrue(baos.toString(), baos.toString().isEmpty());
      assertTrue(mutations.delete());
      assertTrue(!mutant.exists() || FileHelper.deleteAll(mutant));

      //test use mutate count
      assertEquals(errorPrint.toString(), -2, cli.mainInit(new String[] {"-c", "30", "-o", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath(), "-Z"}, NULL_STREAM, errorPrint.printStream()));
      errorPrint.reset();
      assertTrue(baos.toString(), baos.toString().isEmpty());
      assertTrue(!mutations.exists() || mutations.delete());
      assertTrue(!mutant.exists() || FileHelper.deleteAll(mutant));

      //test use mutate rate
      assertEquals(errorPrint.toString(), -3, cli.mainInit(new String[] {"-o", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath(), "-r", "0.5", "-Z"}, NULL_STREAM, errorPrint.printStream()));
      errorPrint.reset();
      assertTrue(baos.toString(), baos.toString().isEmpty());
      assertTrue(!mutations.exists() || mutations.delete());
      assertTrue(!mutant.exists() || FileHelper.deleteAll(mutant));

      //test use twin
      assertEquals(errorPrint.toString(), -4, cli.mainInit(new String[] {"-o", mutant.getPath(), "-O", mutant.getPath(), "-s", mutations.toString(), "-i", in.getPath(), "-Z"}, NULL_STREAM, errorPrint.printStream()));
      errorPrint.reset();
      assertTrue(baos.toString(), baos.toString().isEmpty());

      //      assertTrue(!mutations.exists() || mutations.delete());
      //TODO assertEquals(0, new GenomeMutatorCli().mainInit(new String[] {"--Xmnp-rate", "0.2", "--Xindel-rate", "0.4", "--Xsnp-rate", "0.4", "-c", "1", "-o", mutantDir.getPath(), "-s", mutations.getPath(), "-i",  in.getPath()},
      //     new ByteArrayOutputStream(), NULL_PRINTSTREAM));
      //TODO assertEquals(0, new GenomeMutatorCli().mainInit(new String[] {"--Xmnp-rate", "0.2", "--Xindel-rate", "0.4",
      //          "--Xsnp-rate", "0.4", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(),
      //          "-i",  in.getPath(), "-O", twin.getPath()},
      //          NULL_STREAM, NULL_PRINTSTREAM));

    } finally {
      //final String str = logStream.toString();
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testTabix() throws Exception {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("GenomeMutator", "tmp");
    try {
      final File out = new File(tmp, "out");
      final File in = BuildTestUtils.prereadDNA(mDir, ">a\nacgtacgatcagcatctgac\n");
      final File mutations = new File(tmp, "mutations.txt.gz");
      final GenomeMutatorCli cli = new GenomeMutatorCli();
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int outcode = cli.mainInit(new String[] {"-o", out.getPath(), "-i", in.getPath(), "-s", mutations.getPath(), "--Xmutation-rate", "0.1"}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, outcode);
      assertTrue(new File(mutations.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      TestUtils.containsAll(FileUtils.fileToString(new File(out, "snpsim.log")), "SnpIndex");
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }
}
