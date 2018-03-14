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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SvPatternsCliTest extends AbstractParamsCliTest<BreakpointPatternParams> {
  @Override
  protected ParamsCli<BreakpointPatternParams> getParamsCli() {
    return new SvPatternsCli();
  }

  private static final String EXP_F1 = "Error: You must provide a value for -o DIR" + LS;
  public void testErrorF1() {
    final String err = checkHandleFlagsErr();
    assertTrue("<" + EXP_F1 + "> was not contained in <" + err + ">", err.contains(EXP_F1));
  }

  public void testInitParams() {
    checkHelp("svpatterns [OPTION]... -o DIR FILE+",
        " -o DIR -I FILE ",
        "--region=REGION ", "if set, only process SAM records within the specified range",
        "-I, --input-list-file=FILE file containing a list of VCF format files (1 per line) of discordant breakpoints",
        "--max-fragment-length=INT how far from the breakpoint to look ahead for inversions (Default is 500)",
        "--max-same-distance=INT how far apart can breakpoints be yet still be considered the same place (Default is 50)",
        "--min-support=INT minimum number of supporting reads for a breakend (Default is 3)"
    );
  }
  public void testValidator() throws Exception {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final File in = new File(tmpDir, "discord.vcf.gz");
      final File tmpFile = FileUtils.stringToFile(in.toString(), FileHelper.createTempFile(tmpDir));

      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", in);
      final TabixIndexer tabixIndexer = new TabixIndexer(in, new File(in.getPath() + TabixIndexer.TABIX_EXTENSION));
      tabixIndexer.saveVcfIndex();

      final File outDir = new File(tmpDir, "out");

      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", outDir.getPath(), "Error: No input files specified"));
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-I", tmpFile.getPath(), outDir.getPath(), "Error: You must provide a value for -o DIR"));
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", outDir.getPath(), "-I", tmpFile.getPath(), "--region", "I_am_not_a_region:I_AM_A_FREE_STRING", "is malformed"));

      checkMainInitOk("-o", outDir.getPath(), "-I", tmpFile.getPath());
    }
  }
  public void testMakeParams() throws IOException, UnindexableDataException {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final File in = new File(tmpDir, "discord.vcf.gz");
      final File tmpFile = FileUtils.stringToFile(in.toString(), FileHelper.createTempFile(tmpDir));

      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", in);
      final TabixIndexer tabixIndexer = new TabixIndexer(in, new File(in.getPath() + TabixIndexer.TABIX_EXTENSION));
      tabixIndexer.saveVcfIndex();

      final File outDir = new File(tmpDir, "out");

      final CFlags flags = new CFlags();
      SvPatternsCli.initLocalFlags(flags);

      final String region = "chr1:200-500";
      assertTrue(flags.setFlags("-o", outDir.getPath(), "-I", tmpFile.getPath(), "--region", region, "--max-fragment-length", "90000", "--max-same-distance", "1234", "--min-support", "43"));
      final BreakpointPatternParams params = SvPatternsCli.makeParamsLocal(flags);
      assertEquals(region, params.region().toString());
      assertEquals(1, params.files().size());
      assertEquals(in, params.files().get(0));
      assertEquals(outDir, params.directory());
      assertEquals(90000, params.fragmentLength());
      assertEquals(1234, params.sameDistance());
      assertEquals(43, params.minDepth());
    }
  }
  public void testMakeParamsDefaults() throws IOException, UnindexableDataException {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final File in = new File(tmpDir, "discord.vcf.gz");
      final File tmpFile = FileUtils.stringToFile(in.toString(), FileHelper.createTempFile(tmpDir));

      FileHelper.resourceToFile("com/rtg/variant/sv/discord/pattern/resources/discordant_pairs.vcf.gz", in);
      final TabixIndexer tabixIndexer = new TabixIndexer(in, new File(in.getPath() + TabixIndexer.TABIX_EXTENSION));
      tabixIndexer.saveVcfIndex();

      final File outDir = new File(tmpDir, "out");

      final CFlags flags = new CFlags();
      SvPatternsCli.initLocalFlags(flags);

      assertTrue(flags.setFlags("-o", outDir.getPath(), "-I", tmpFile.getPath()));
      final BreakpointPatternParams params = SvPatternsCli.makeParamsLocal(flags);
      assertEquals(null, params.region());
      assertEquals(500, params.fragmentLength());
      assertEquals(50, params.sameDistance());
      assertEquals(3, params.minDepth());
    }
  }
}
