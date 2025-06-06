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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
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
    for (int i = 0; i < 100; ++i) {
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
        , "-t,", "--template=SDF", "reference genome"
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

        checkHandleFlagsErr();
        checkHandleFlagsErr("-i", sam.getAbsolutePath());
        checkHandleFlagsErr("-t", template.getAbsolutePath());
        checkHandleFlagsErr("-t", template.getAbsolutePath(), "-I", "blah");

        checkMainInitOk("-t", template.getAbsolutePath(), sam.getAbsolutePath());
        checkMainInitOk("-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", template.getAbsolutePath());
        checkMainInitOk("-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", template.getAbsolutePath(), "--Xignore-cg-fragment");

        MainResult res = MainResult.run(getCli(), "-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", pairedReads.getAbsolutePath());
        assertEquals(res.err(), 0, res.rc());
        assertFalse(res.out().contains("Distribution"));

        res = MainResult.run(getCli(), "-t", template.getAbsolutePath(), sam.getAbsolutePath(), "-r", pairedReads.getAbsolutePath(), "-D");
        assertEquals(res.err(), 0, res.rc());
        assertTrue(res.out().contains("Distribution of read hits"));
      } finally {
        assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
      }
    } finally {
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
    checkMainInitOk("-D", "--validate", "-t", template.getPath(), "--Xmismatch-penalty", "1", "--Xgap-open-penalty", "1", "--Xgap-extend-penalty", "1", "--Xunknowns-penalty", "0", sam.getPath());
  }

}
