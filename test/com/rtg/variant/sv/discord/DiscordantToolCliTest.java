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
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class DiscordantToolCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new DiscordantToolCli();
  }

  public void testFlagValidator() throws IOException {
    String str = checkHandleFlagsErr("-t", "genome-stupid-name", "-o", "output", "-r", "blah.txt", "blah.txt");
    assertTrue(str, str.contains("The specified SDF, \"genome-stupid-name\", does not exist."));
    final File tempDir = FileUtils.createTempDir("discordantTool", "flagValidatorTest");
    try {
      str = checkHandleFlagsErr("-t", tempDir.getPath(), "-o", "randomdir_output", "-r", "blah.txt", "blah.txt", "-s", "0");
      assertTrue(str, str.contains("Expected a positive integer for parameter \"min-support\"."));
      str = checkHandleFlagsOut("-t", tempDir.getPath(), "-o", "randomdir_output", "-r", "blah.txt", "blah.txt");
      assertEquals("", str);
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testMakeParams() throws InvalidParamsException, IOException {
    final File tempDir = FileUtils.createTempDir("discordantTool", "makeParamsTest");
    try {
      final File f = new File(tempDir, "blah.txt");
      FileUtils.stringToFile("RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", f);
      final File gen = ReaderTestUtils.getDNADir(">g1" + LS + "aaatcgactggtcagctagg" + LS, tempDir);
      final CFlags flags = new CFlags();
      DiscordantToolCli.initFlags(flags);
      flags.setName("blah");
      flags.setInvalidFlagHandler(CFlags.DEFAULT_INVALID_FLAG_HANDLER);
      flags.setFlags("--template", gen.getPath(), "--output", new File(tempDir, "blah").getPath(), "-r", f.getPath(), f.getPath());
      final DiscordantToolParams params = DiscordantToolCli.makeParams(flags);
      assertTrue(params.genome().reader() instanceof DefaultSequencesReader);
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
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
//        , "--max-coverage=INT", "if set, will only output variants where coverage is less than this amount"
//        , "--max-ambiguity=INT", "threshold for ambiguity above which calls are not made"
        );
    checkExtendedHelp("rtg discord"
        , "--Xdebug-output", "produce debug output in addition to VCF"
        );
  }
}
