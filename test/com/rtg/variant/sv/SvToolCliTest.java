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
package com.rtg.variant.sv;

import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.SAM_UNSORTED;
import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SvToolCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SvToolCli();
  }

  /** Test for an error that will be picked up during flags parsing. */
  private void checkFlagsError(final String[] args, final String exp) {
    assertTrue(checkHandleFlagsErr(args).contains(exp));
  }

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -t SDF" + LS;

  public void testErrorFlags() {
    checkFlagsError(new String[] {}, EXP_F1);
  }

  private void checkParamsError(final String[] args, final String... exp) {
    TestUtils.containsAll(checkHandleFlagsErr(args).replaceAll(LS, " ").replaceAll("\\s+", " "), exp);
  }

  private static final String EXP_P1 = "The specified SDF, \"blah\", does not exist";
  private static final String[] EXP_P3 = {"The specified file, ", ", is not an SDF"};

  public void testErrorParams() throws IOException, InvalidParamsException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final File out = new File(tempDir, "out");
      checkParamsError(new String[]{"-o", out.getPath(), "-t", "blah", "blah.alignments"}, EXP_P1);
      FileHelper.deleteAll(out);
      final File tmp = File.createTempFile("test", "tmp", tempDir);
      final File aln = new File(tempDir, "blah.alignments");
      checkParamsError(new String[]{"-o", out.getPath(), "-t", tmp.getPath(), aln.getPath()}, EXP_P3);
      final File gen = FileUtils.createTempDir("test", "genome", tempDir);
      assertTrue(aln.createNewFile());
      final File rg = new File(tempDir, "rgfile.txt");
      assertTrue(rg.createNewFile());
      checkHandleFlagsOut("-o", out.getPath(), "-t", gen.getPath(), aln.getPath(), "-r", rg.getPath());

      assertTrue(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-b", "0").contains("Expected a positive integer for parameter \"Xbin-size\""));
      assertTrue(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-s", "0").contains("Expected a positive integer for parameter \"step\""));
      assertTrue(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-f", "0").contains("Expected a positive integer for parameter \"fine-step\""));
      assertTrue(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-s", "1").contains("Parameter \"fine-step\" should be smaller than or equal to parameter \"step\""));
    }
  }

  public void testMakeParams() throws IOException, InvalidParamsException {
    SvToolTaskTest.check(REF_SEQS, SAM_UNSORTED, "RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new String[] {}, null, new String[] {"is not sorted in coordinate order."}, "", 1, true, null);
  }

  public void testInitParams() {
    checkHelp("rtg sv"
        , "Analyses SAM records to determine the location of structural variants."
        , "[OPTION]... -o DIR -t SDF -I FILE -R FILE"
        , "[OPTION]... -o DIR -t SDF -r FILE FILE+"
        , "simple-signals", "if set, also output simple signals"
        , "-f,", "--fine-step=INT", "step size in interesting regions (Default is 10)"
        , "-s,", "--step=INT", "step size (Default is 100)"
        );
    checkExtendedHelp("rtg sv"
        , "-b,", "--Xbin-size=INT", "bin size used by simple signals (Default is 10)"
        , "--Xcorrections=FILE", "file containing per position corrections"
        , "--Xheterozygous", "if set, also include heterozygous bayesian models"
        );
  }
}
