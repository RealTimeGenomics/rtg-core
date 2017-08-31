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

package com.rtg.variant.sv.discord;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.DefaultSequencesReader;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 */
public class DiscordantToolCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new DiscordantToolCli();
  }

  public void testFlagValidator() throws IOException {
    TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", "genome-name", "-o", "output", "-r", "blah.txt", "blah.txt"),
      "The specified SDF, \"genome-name\", does not exist.");
    try (final TestDirectory tempDir = new TestDirectory()) {
      final File f = new File(tempDir, "blah.txt");
      assertTrue(f.createNewFile());
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", tempDir.getPath(), "-o", "randomdir_output", "-r", f.getPath(), f.getPath(), "-s", "0"),
        "--min-support", "must be at least 1");
      checkHandleFlagsOut("-t", tempDir.getPath(), "-o", "randomdir_output", "-r",  f.getPath(), f.getPath());
    }
  }

  public void testMakeParams() throws InvalidParamsException, IOException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final File f = new File(tempDir, "blah.txt");
      assertTrue(f.createNewFile());
      FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", f);
      final File gen = ReaderTestUtils.getDNADir(">g1" + LS + "aaatcgactggtcagctagg" + LS, tempDir);
      checkHandleFlags("--template", gen.getPath(), "--output", new File(tempDir, "blah").getPath(), "-r", f.getPath(), f.getPath(), "--min-support=17", "--overlap-fraction=2.3");
      final DiscordantToolParams params = ((DiscordantToolCli) mCli).makeParams();
      assertTrue(params.genome().reader() instanceof DefaultSequencesReader);
      assertEquals(17, params.minBreakpointDepth());
      assertEquals(2.3, params.overlapFraction());
    }
  }

  public void testInitFlags() {
    checkHelp("rtg discord"
        , "Analyses SAM records to determine the location of breakends."
        , "[OPTION]... -o DIR -t SDF -I FILE -R FILE"
        , "[OPTION]... -o DIR -t SDF -r FILE FILE+"
        , "--bed", "produce output in BED format in addition to VCF"
        , "--consistent-only", "only include breakends with internally consistent supporting reads"
        , "-s,", "--min-support=INT", "minimum number of supporting reads for a breakend (Default is 3)"
        , "-c,", "--max-hits=INT", "if set, ignore SAM records with an alignment count that exceeds this value"
        );
    checkExtendedHelp("rtg discord", "--Xdebug-output", "produce debug output in addition to VCF", "--Xmultisample");
  }
}
