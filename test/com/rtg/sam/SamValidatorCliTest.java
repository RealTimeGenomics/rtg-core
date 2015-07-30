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
package com.rtg.sam;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 * Tests corresponding class.
 */
public class SamValidatorCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SamValidatorCli();
  }

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

  private String getSequence() {
    final StringBuilder sb = new StringBuilder(200);
    sb.append(">a\n");
    for (int i = 0; i < 100; i++) {
      sb.append("acgtacgatcagcatctgac");
    }
    sb.append("\n");
    return sb.toString();
  }

  public void testHelp() {
    checkHelp("rtg samstats"
        , "Prints alignment statistics from the contents of the output SAM/BAM file."
        , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
        , "-r,", "--reads=SDF", "reads SDF"
        , "-t,", "--template=SDF", "template SDF"
        , "FILE+", "SAM/BAM result file (must contain read-ids not read names). May be specified 0 or more times"
        , "--consensus", "record consensus data. Requires roughly 5 times template length of RAM"
        , "-D,", "--distributions", "display distributions of insert sizes, alignment scores and read hits"
        , "--per-file", "output per-file statistics"
        , "--validate", "validate mapping of read to reference. Tests matching of bases according to CIGAR format"
        );
    checkExtendedHelp("rtg samstats"
        , "--Xignore-cg-fragment", "ignore unusual Complete Genomics fragment lengths."
        );
  }

  public void testFlagValidator() throws IOException {
    final File template = ReaderTestUtils.getDNASubDir(getSequence(), mDir);
    final File pairedReads = FileUtils.createTempDir("temp", "paired", mDir);
    final SdfId sdfId = ReaderTestUtils.createPairedReaderDNA(getSequence(), getSequence(), pairedReads, null);
    final File sam = new File(mDir, "samvalnosort.sam");
    final String tab = "\t";
    final String samHeader = ""
      + "@HD" + tab + "VN:1.0" + tab + "SO:coordinate" + StringUtils.LS
      + "@SQ" + tab + "SN:a" + tab + "LN:30" + StringUtils.LS
      + "@CO" + tab + "READ-SDF-ID:" + sdfId + StringUtils.LS
      + "0" + tab + "0" + tab + "a" + tab + "2" + tab + "255" + tab + "10M" + tab + "*" + tab + "0" + tab + "0" + tab + "AAAAAAAAAA" + tab +  "IB7?*III<I" + tab + "AS:i:0" + tab + "IH:i:1" + StringUtils.LS
      ;

    FileUtils.stringToFile(samHeader, sam);
    final File inputList = new File(mDir, "input.txt");
    assertTrue(inputList.createNewFile());
    checkHandleFlagsErr("-t", template.getPath(), "-I", inputList.getPath());
    assertTrue(inputList.delete());
    FileUtils.stringToFile(sam.getPath() + StringUtils.LS, inputList);
    String err = checkHandleFlagsErr("-t", "blah", "-I", inputList.getPath());
    assertTrue(err, err.contains("The specified SDF, \"blah\", does not exist."));
    err = checkHandleFlagsErr("-t", sam.getPath(), "-I", inputList.getPath());
    assertTrue(err, err.contains("The specified file, \"" + sam.getPath() + "\", is not an SDF."));
    err = checkHandleFlagsErr("-t", template.getPath(), "-I", inputList.getPath(), "-r", "blah");
    assertTrue(err, err.contains("The specified SDF, \"blah\", does not exist."));
    err = checkHandleFlagsErr("-t", template.getPath(), "-I", inputList.getPath(), "-r", sam.getPath());
    assertTrue(err, err.contains("The specified file, \"" + sam.getPath() + "\", is not an SDF."));
    err = checkHandleFlagsErr("-t", template.getPath(), "-I", inputList.getPath(), "-r", pairedReads.getPath(), "--Xmismatch-penalty", "9");
    assertTrue(err, err.contains("Must specify all of --Xmismatch-penalty, --Xgap-open-penalty, --Xgap-extend-penalty, --Xunknowns-penalty if specifying any."));
    checkHandleFlagsOut("-t", template.getPath(), "-I", inputList.getPath(), "-r", pairedReads.getPath());
  }

  private static final String TAB = "\t";
  private static final String SAM_HEAD = ""
        + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + StringUtils.LS
        + "@SQ" + TAB + "SN:a" + TAB + "LN:30" + StringUtils.LS;

  public void testFlags() throws IOException {
    ByteArrayOutputStream outbaos = new ByteArrayOutputStream();
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream errStr = new PrintStream(baos);
    final File template = ReaderTestUtils.getDNASubDir(getSequence(), mDir);
    final File pairedReads = FileUtils.createTempDir("temp", "paired");
    final SdfId sdfId = ReaderTestUtils.createPairedReaderDNA(getSequence(), getSequence(), pairedReads, null);

    try {
      final File dir = FileUtils.createTempDir("samvaltest", "samfileandrecord");
      final File sam = new File(dir, "samvalnosort.sam");
      try {
        final String samHeader = ""
          + SAM_HEAD
          + "@CO" + TAB + "READ-SDF-ID:" + sdfId + StringUtils.LS
          + "0" + TAB + "0" + TAB + "a" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + StringUtils.LS
          ;

        FileUtils.stringToFile(samHeader, sam);

        int result;
        result = new SamValidatorCli().mainInit(new String[] {}, outbaos, errStr);
        assertEquals(1, result);
        result = new SamValidatorCli().mainInit(new String[] {"-h"}, outbaos, errStr);
        outbaos.flush();
        assertTrue(outbaos.toString().contains("Usage:"));
        assertEquals(1, result);
        result = new SamValidatorCli().mainInit(new String[] {"-i", sam.getAbsolutePath()}, outbaos, errStr);
        assertEquals(1, result);
        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath()}, outbaos, errStr);
        assertEquals(1, result);
        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), "-I", "blah"}, outbaos, errStr);
        assertEquals(1, result);

        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), sam.getAbsolutePath()}, outbaos, errStr);
        assertEquals(0, result);

        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", template.getAbsolutePath()}, outbaos, errStr);
        assertEquals(0, result);

        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", template.getAbsolutePath(), "--Xignore-cg-fragment"}, outbaos, errStr);
        assertEquals(0, result);

        outbaos = new ByteArrayOutputStream();
        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", pairedReads.getAbsolutePath()}, outbaos, errStr);
        assertEquals(0, result);
        outbaos.flush();
        assertFalse(outbaos.toString().contains("Distribution"));

        outbaos = new ByteArrayOutputStream();
        result = new SamValidatorCli().mainInit(new String[] {"-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", pairedReads.getAbsolutePath(), "-D"}, outbaos, errStr);
        assertEquals(0, result);
        outbaos.flush();
        final String out = outbaos.toString();
        assertTrue(out.contains("Distribution of read hits"));
      } finally {
        assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
      }
    } finally {
      assertTrue(!template.exists() || FileHelper.deleteAll(template));
      assertTrue(!pairedReads.exists() || FileHelper.deleteAll(pairedReads));
    }
  }

  public void testEndToEndCG() throws Exception {
    final File template = ReaderTestUtils.getDNASubDir(">t\nGGATTGAGACTGGTAAAATATGAAGTGACCACCAAAGGGAGCTTGAGAGA\n", mDir);

    final String samstr = ""
            + SAM_HEAD
            + "0" + TAB + "179" + TAB + "t" + TAB + "5" + TAB + "55" + TAB + "1=1X4=1X17=6N10=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TNAGACNGGTAAAATATGAAGTGAAAGGGAGCTT" + TAB +  "1!551-!//2-+2,/22002535./413664212" + TAB + "AS:i:31" + TAB + "IH:i:1" + TAB + "XU:Z:1=1R3=1B1X1=1R17=6N10=" + TAB + "XQ:Z:+" + TAB + "XR:Z:T" + StringUtils.LS
            ;

    final File sam = FileHelper.stringToGzFile(samstr, new File(mDir, "blah.sam"));

    final SamValidatorCli svc = new SamValidatorCli();

    final MemoryPrintStream omps = new MemoryPrintStream();
    final MemoryPrintStream emps = new MemoryPrintStream();

    final int code = svc.mainInit(new String[] {"-D", "--validate", "-t", template.getPath(), "--Xmismatch-penalty", "1", "--Xgap-open-penalty", "1", "--Xgap-extend-penalty", "1", "--Xunknowns-penalty", "0", sam.getPath()}, omps.outputStream(), emps.printStream());
    assertEquals(emps.toString(), 0, code);

    assertEquals(emps.toString(), 0, emps.toString().length());
  }

}
