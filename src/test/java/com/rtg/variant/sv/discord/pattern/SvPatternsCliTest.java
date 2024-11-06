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
