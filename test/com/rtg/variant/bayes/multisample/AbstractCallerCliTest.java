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

package com.rtg.variant.bayes.multisample;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.Environment;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantParams;

/**
 */
public abstract class AbstractCallerCliTest extends AbstractCliTest {

  protected abstract String getModuleName();

  private String spaces(final int indent) {
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < indent; k++) {
      sb.append(' ');
    }
    return sb.toString();
  }

  /** Test for an error that will be picked up during params object construction. */
  public void checkParamsError(final String[] args0, ErrorType expErrorType) throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    final File subjectsDir = FileHelper.createTempDirectory();
    try {
      ReaderTestUtils.getReaderDNA(">t\nacgt", subjectsDir, null).close();
      final File output = FileUtils.createTempDir("fileout", "varianceCli");
      try {
        final File map = new File(output, "map.sam.gz");
        FileHelper.stringToGzFile(SharedSamConstants.SAM9, map);
        assertTrue(new File(map.getPath() + TabixIndexer.TABIX_EXTENSION).createNewFile());
        final String sub = subjectsDir.getPath();
        final String[] args = Utils.append(args0, "-t", sub, "-o", output.getPath(), map.getPath(), "-m", "default");
        final AbstractMultisampleCli v = (AbstractMultisampleCli) mCli;
        checkHandleFlags(args);
        try {
          v.makeParams();
          fail();
        } catch (final InvalidParamsException e) {
          assertEquals(expErrorType, e.getErrorType());
//          assertTrue("Exception:" + e.getClass().getName() + " Actual: " + e.getMessage() + "\n" + "Expected to contain: " + exp, e.getMessage().contains(exp));
        }
      } finally {
        assertTrue(FileHelper.deleteAll(output));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }
  }

  /** Test for an error that will be picked up during params object construction. */
  public void checkFlagVaidations(final String[] args0, final String exp) throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    final File subjectsDir = FileHelper.createTempDirectory();
    try {
      ReaderTestUtils.getReaderDNA(">t\nacgt", subjectsDir, null).close();
      final File output = FileUtils.createTempDir("fileout", "varianceCli");
      try {
        final File map = new File(output, "map.sam.gz");
        FileHelper.stringToGzFile(SharedSamConstants.SAM9, map);
        assertTrue(new File(map.getPath() + TabixIndexer.TABIX_EXTENSION).createNewFile());
        final String sub = subjectsDir.getPath();
        final String[] args = Utils.append(args0, "-t", sub, "-o", new File(output, "snpscalls").getPath(), map.getPath(), "-m", "default");
        final String res = checkHandleFlagsErr(args);
        TestUtils.containsAll(res, exp);
      } finally {
        assertTrue(FileHelper.deleteAll(output));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }
  }

  public void testInitParams() {
    checkHelp(getModuleName() + " [OPTION]... -o DIR -t SDF FILE+",
        "SAM/BAM format files containing mapped reads",
        "directory for output",
        "of the reference genome the reads",
        "write variant calls covering every position irrespective of thresholds",
        "do not gzip the output",
        "for mated reads that have no mapping quality supplied",
        "for unmated reads that have no mapping quality supplied",
        "if set, ignore mapping calibration files",
        "SDF of the reference genome the reads have been mapped against",
        "if set, force sequencer machine settings. One of [default, illumina, ls454_se, ls454_pe, complete, iontorrent]",
        "if set, will output simple SNPs only",
        "for mated reads that have no mapping quality supplied use this as the default quality (in Phred format from 0 to 63)",
        "for unmated reads that have no mapping quality supplied use this as the default quality (in Phred format from 0 to 63)",
        "Calls sequence variants",
        "file containing a list of SAM/BAM format files (1 per line) containing mapped reads",
        "ploidy to use",
        "sex of individual"
        );
  }

  private final String mUsagePrefix = ""
      + "Usage: rtg " + getModuleName() + " [OPTION]... -o DIR -t SDF FILE+" + LS
      + "           " + spaces(getModuleName().length()) + " [OPTION]... -o DIR -t SDF -I FILE" + LS;

  private static final String USAGE_SUFFIX = ""
      + "" + LS
      + "Try '--help' for more information"
      ;

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -t SDF" + LS;

  public void testErrorF1() throws InvalidParamsException {
    assertTrue(checkHandleFlagsErr().contains(EXP_F1));
  }

  //Error in prior file
  public void testErrorP1() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--Xpriors", "testpriorbad1"}, ErrorType.PROPS_KEY_NOT_FOUND);
  }

  //Error in prior file
  public void testErrorP2() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--Xpriors", "testpriorbad2"}, ErrorType.PRIOR_KEY_VALUE_INVALID);
  }

  //Error in prior file
  public void testErrorP3() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--Xpriors", "testpriorbad3"}, ErrorType.PRIOR_KEY_VALUE_INVALID);
  }

  //Error in prior file
  public void testErrorP4() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--Xpriors", "testpriorbad4"}, ErrorType.PROPS_INVALID);
  }

  //Prior file missing
  public void testErrorP5() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--Xpriors", "testpriordoesntexist"}, ErrorType.INFO_ERROR);
  }

  //Prior file missing
  public void testErrorP6() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"-q", "64"}, ErrorType.INVALID_INTEGER);
  }

  public void testMaxCoverageAndCoverageMultiplierSet() throws InvalidParamsException, IOException {
    checkFlagVaidations(new String[] {"--filter-depth", "64", "--filter-depth-multiplier", "4.0"}, "Only one of --filter-depth or --filter-depth-multiplier can be set");
  }

  //Prior file missing
  public void testErrorP7() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"-q", "-1"}, ErrorType.INVALID_INTEGER);
  }

  //Prior file missing
  public void testErrorP8() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"--rdefault-mated", "64"}, ErrorType.INVALID_INTEGER);
  }

  private static final String SAM = "@HD\tVN:1.0\tSO:coordinate\n"
      + "@SQ\tSN:t0\tLN:4\n"
      + "0\t16\t*\t0\t0\t*\t*\t0\t0\tAGCTTGTCATTCTGACTGCAACGGGCAATATG\t*";

  private VariantParams getParams(final String[] args0) throws InvalidParamsException, IOException {
    final File subjectsDir = FileUtils.createTempDir("test", "variant");
    ReaderTestUtils.getReaderDNA(">t0\nacgt", subjectsDir, null).close();
    final String sub = subjectsDir.getPath();
    final String[] args = Utils.append(args0, "-t", sub);
    checkHandleFlagsOut(args);
    final VariantParams vp = ((AbstractMultisampleCli) mCli).makeParams();
    vp.integrity();
    vp.close();
    return vp;
  }

  public void testNotExistingInput() throws Exception {
    Diagnostic.setLogStream();
    final File genome = ReaderTestUtils.getDNADir(">genome\nacgt");
    try {
      final ByteArrayOutputStream ba = new ByteArrayOutputStream();
      final PrintStream err = new PrintStream(ba);
      final File outFile = FileUtils.createTempDir("reallyreallydoesntexist", "tempout");
      assertTrue(FileHelper.deleteAll(outFile));
      try {
        final String[] args = {"-o", outFile.getPath(), "-t", genome.getPath(), "output" + StringUtils.FS + "sam.gz"};
        final int res = getCli().mainInit(args, new ByteArrayOutputStream(), err);
        assertEquals(1, res);
        err.close();
        TestUtils.containsAll(ba.toString(), "Error: File not found: \"output" + StringUtils.FS + "sam.gz\"",
          "Error: There were 1 invalid input file paths");
      } finally {
        assertTrue(FileHelper.deleteAll(outFile));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(genome));
    }
  }

  public void testNotExistingGenomeDir() throws Exception {
    Diagnostic.setLogStream();
    final File output = FileUtils.createTempDir("carianceclitest", "nonexistinggenome");
    FileHelper.deleteAll(output);
    try {
      final ByteArrayOutputStream ba = new ByteArrayOutputStream();
      final PrintStream err = new PrintStream(ba);
      final String[] args = {"-t", "input", "-o", output.getPath(), new File(output, "sam.gz").getPath()};
      final int res = getCli().mainInit(args, new ByteArrayOutputStream(), err);
      assertEquals(1, res);
      err.close();
      final String ex = "Error: The specified SDF, \"input\", does not exist." + LS + mUsagePrefix + USAGE_SUFFIX;
      // assertEquals(ex, flags.getInvalidFlagMsg());
      assertEquals(ex + LS, ba.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(output));
    }
  }

  private static final String TEMPALTE = ">chr1" + LS
      + "aaaaaaaaccccccccgggggggg";


  private static final String SAM_INPUT2 = ""
      + "@HD\tVN:1.0\tSO:coordinate" + LS
      + "@PG\tID:rtg\tVN:vPOST-2.0.1-DEV build 31834 (2010-09-30)\tCL:map" + LS
      + "@SQ\tSN:chr1\tLN:24" + LS
      ;

  public void testSafeSexWarning() throws Exception {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("carianceclitest", "nonexistinggenome");
    try {
      final File template = new File(f, "template");
      ReaderTestUtils.getReaderDNA(TEMPALTE, template, null).close();
      final File out = new File(f, "out");
      final File sam = new File(f, "x.sam");
      FileUtils.stringToFile(SAM_INPUT2, sam);
      try (MemoryPrintStream err = new MemoryPrintStream()) {
        final String[] args = {"-t", template.getPath(), "-o", out.getPath(), sam.getPath(), "--sex", "male"};
        final int res = getCli().mainInit(args, new ByteArrayOutputStream(), err.printStream());
        assertEquals(1, res);
        final String e = err.toString();
        assertTrue(e, e.startsWith("Sex-specific processing was specified but"));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  protected VariantParams params(final String...args) throws InvalidParamsException, IOException, UnindexableDataException {
    Diagnostic.setLogStream();
    final AbstractMultisampleCli vc = (AbstractMultisampleCli) getCli();
    final File outFile = FileUtils.createTempDir("variancetest", "tempout");
    FileHelper.deleteAll(outFile);
    final File inFile = FileUtils.createTempDir("variancetest", "tempin");
    try {
      final File map = new File(inFile, "map.gz");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(SAM.getBytes()), map);
      new TabixIndexer(map, new File(inFile, "map.gz.tbi")).saveSamIndex();

      final String[] fullArgs = com.rtg.util.Utils.append(new String[]{"-o", outFile.getPath(), map.getPath(), "-m", "default"}, args);
      final VariantParams vp = getParams(fullArgs);
      try {
        final OutputStream out = new ByteArrayOutputStream();

        final IORunnable task = vc.task(vp, out);
        assertNotNull(task);
        task.run();
        return vp;
      } finally {
        vp.close();
        FileHelper.deleteAll(vp.genome().directory());
        FileHelper.deleteAll(vp.directory());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(inFile));
    }
  }

  public void testNonIdentity() throws InvalidParamsException, IOException, UnindexableDataException {
    final VariantParams vp = params();
    assertTrue(vp.nonidentityPosterior());
  }

  public void testDefault() throws InvalidParamsException, IOException, UnindexableDataException {
    final VariantParams vp = params();
    assertEquals(Environment.getAvailableProcessors(), vp.execThreads());
    assertEquals(Environment.getAvailableProcessors(), vp.ioThreads());
    assertEquals(Integer.MAX_VALUE, vp.maxCoverageFilter().thresholdSingle("some sequence"));
    assertTrue(vp.nonidentityPosterior());
  }

  public void testThreads() throws InvalidParamsException, IOException, UnindexableDataException {
    final VariantParams params = params("-T", "13", "--Xio-threads", "7");
    assertEquals(13, params.execThreads());
    assertEquals(7, params.ioThreads());
  }
}
