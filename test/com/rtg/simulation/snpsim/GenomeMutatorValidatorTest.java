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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;


/**
 */
public class GenomeMutatorValidatorTest extends AbstractCliTest {

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
    return new GenomeMutatorCli();
  }

  public void testCliValidator1() throws IOException, InvalidParamsException {
    final File tempDir = FileUtils.createTempDir("genomemutatorclitest", "checkparamserror");
    try {
      final File genomeDir = FileHelper.createTempDirectory();
      try {
        ReaderTestUtils.getReaderDNA("", genomeDir, null).close();
        final File xgenomeFile = new File(tempDir, "xgenome");
        final File mutGenome = new File(tempDir, "mut1");
        final File snps = new File(tempDir, "snps");
        final File snpsZip = new File(tempDir, "snps.gz");
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath()).contains("You must provide values for -o SDF -s FILE"));
        assertTrue(checkMainInitBadFlags("-i", xgenomeFile.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath()).contains("The specified SDF, \"" + xgenomeFile.getPath() + "\", does not exist."));
        FileUtils.stringToFile("kgajkgf", xgenomeFile);
        assertTrue(checkMainInitBadFlags("-i", xgenomeFile.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath()).contains("The specified file, \"" + xgenomeFile.getPath() + "\", is not an SDF."));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "-p", "human", "--Xmutation-count=100").contains("Cannot set --Xpriors and --Xmutation-count"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "-p", "human", "--Xmutation-rate=1").contains("Cannot set --Xpriors and --Xmutation-rate"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xmutation-count=100", "--Xmutation-rate=1").contains("Must specify at most one of --Xmutation-count and --Xmutation-rate"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xmutation-count=-1").contains("Mutation count must be positive"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xmutation-rate=2").contains("Mutation rate must be between 0.0 and 1.0"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xsnp-rate=-1").contains("SNP rate is under 0"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xmnp-rate=-1").contains("MNP rate is under 0"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "--Xindel-rate=-1").contains("INDEL rate is under 0"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "-p", "human", "--Xsnp-rate=1").contains("Cannot set --Xpriors and --Xsnp-rate"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "-p", "human", "--Xindel-rate=1").contains("Cannot set --Xpriors and --Xindel-rate"));
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(),
            "-p", "human", "--Xmnp-rate=1").contains("Cannot set --Xpriors and --Xmnp-rate"));
        FileUtils.stringToFile("pseudo snps", snps);
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath(), "-Z").contains("SNP mapping file already exists"));
        FileUtils.stringToFile("pseudo snps", snpsZip);
        assertTrue(checkMainInitBadFlags("-i", genomeDir.getPath(), "-s", snps.getPath(), "-o", mutGenome.getPath()).contains("SNP mapping file already exists"));
      } finally {
        assertTrue(FileHelper.deleteAll(genomeDir));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  // copied from variant CLI tests
  public void testNotExistingInput() throws Exception {
    final File tmp = FileUtils.createTempDir("GenomeMutator", "tmp");
    try {
      Diagnostic.setLogStream();
      final File outFile = FileUtils.createTempDir("reallyreallydoesntexist", "tempout");
      try {
        assertTrue(FileHelper.deleteAll(outFile));
        final File snpFile = File.createTempFile("genomemutatorcli", "snpout", tmp);
        assertTrue(checkMainInitBadFlags("-i", "missing", "-o", outFile.getPath(), "-s", snpFile.getPath()).contains("The specified SDF, \"missing\", does not exist."));
      } finally {
        FileHelper.deleteAll(outFile);
      }
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testFlags() throws IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);

    final File tmp = FileUtils.createTempDir("GenomeMutator", "tmp");
    try {
      final File mutant = File.createTempFile("GenomeMutator", "mutant", tmp);
      final File mutations = File.createTempFile("GenomeMutator", "mutations", tmp);
      final File twinmutations = File.createTempFile("GenomeMutator", "flagsMutationsTwin", tmp);
      final File twin = File.createTempFile("GenomeMutator", "twin", tmp);
      assertTrue(mutant.delete());
      assertTrue(twin.delete());
      assertTrue(twinmutations.delete());
      final File in = ReaderTestUtils.getDNADir(mDir);
      checkHandleFlagsErr();
      checkHandleFlagsErr(in.getPath());
      checkHandleFlagsErr("-r", "0.5", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i", in.getPath());
      checkHandleFlagsErr("-r", "1.2", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i", in.getPath());
      checkHandleFlagsErr("-r", "-1.2", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i", in.getPath());
      checkHandleFlagsErr("-c", "-1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i", in.getPath());
      checkHandleFlagsErr("--Xindel-rate", "-1", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i", in.getPath());
      checkHandleFlagsErr("--Xsnp-rate", "-1", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xmnp-rate", "-1", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xindel-rate", "1.2", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xsnp-rate", "1.2", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xmnp-rate", "1.2", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xsnp-rate", "0.6", "--Xmnp-rate", "0.6", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xindel-rate", "0.6", "--Xmnp-rate", "0.6", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      checkHandleFlagsErr("--Xindel-rate", "0.6", "--Xsnp-rate", "0.6", "-c", "1", "-o", mutant.getPath(), "-s", mutations.getPath(), "-i",  in.getPath());
      assertTrue(!mutations.exists() || mutations.delete());
    } finally {
      //final String str = logStream.toString();
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }
}
