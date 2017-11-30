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
package com.rtg;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class CoreCommandTest extends TestCase {

  @Override
  public void setUp() {
    GlobalFlags.resetAccessedStatus();
  }

  public void testUsage() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    ToolsCommand.FORMAT.mainInit(new String[]{"-h"}, baos, null);
    baos.flush();
    TestUtils.containsAllUnwrapped(baos.toString(), "Usage: rtg format [OPTION]... -o SDF FILE+",
        "-h, --help", "print help on command-line flag usage");
    assertEquals("FORMAT", ToolsCommand.FORMAT.getCommandName());
    assertEquals("FORMAT", ToolsCommand.FORMAT.toString());
    assertEquals(CommandCategory.FORMAT, ToolsCommand.FORMAT.getCategory());
  }

  public void testModules() {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(baos);
    assertEquals(0, ToolsCommand.VERSION.mainInit(new String[0], baos, ps));
    assertEquals(0, CoreCommand.LICENSE.mainInit(new String[0], baos, ps));
    assertEquals(0, CoreCommand.HELP.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CG2SDF.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CGMAP.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.CGSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNV.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CALIBRATE.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.MENDELIAN.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.COVERAGE.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CHRSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.FORMAT.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.GENOMESIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPX.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.COMPOSITIONMETA.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.FUNCTIONALMETA.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.METAGENOMICS.mainInit(new String[0], baos, ps));

    assertEquals(1, CoreCommand.SIMILARITY.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAP.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPF.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.READSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAM2BAM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMRENAME.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMMERGE.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPXRENAME.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SDF2FASTA.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SDF2FASTQ.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDF2QUALA.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SDFSUBSET.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SDFSUBSEQ.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDFSPLIT.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SDFSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFFILTER.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFANNOTATE.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SNPINTERSECT.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.TAXFILTER.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.NCBI2TAX.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.TAXSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFMERGE.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.POPSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SAMPLESIM.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.CHILDSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.DENOVOSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.SAMPLEREPLAY.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFDECOMPOSE.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFEVAL.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SINGLETON.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.READSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSIM.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SEGMENT.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSUMMARY.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVPANELBUILD.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.READSIMEVAL.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.ASSEMBLE.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.ADDPACBIO.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.INDEX.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVIEW.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SPECIES.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.METASNP.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SVPREP.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SV.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.DISCORD.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_FAMILY.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_SOMATIC.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_POPULATION.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.BGZIP.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.EXTRACT.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRBUILD.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRPREDICT.mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.PEDFILTER.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.PEDSTATS.mainInit(new String[0], baos, ps));
    assertEquals(1, ToolsCommand.VCFSUBSET.mainInit(new String[0], baos, ps));
  }

  public void testModuleHelp() {
    final String[] helpArgs = {"--help"};
    for (final Command module : CoreCommand.INFO.commands()) {
      if ((module == ToolsCommand.VERSION) || (module == CoreCommand.LICENSE) || (module == CoreCommand.HELP)) {
        continue; // These modules have no help
      }
      GlobalFlags.resetAccessedStatus();
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      final PrintStream ps = new PrintStream(baos);
      module.mainInit(helpArgs, baos, ps);
      assertTrue(module.getCommandName() + " produces no help output: " + baos.toString(), baos.toString().contains("Usage"));
    }
  }

  public void testReleaseLevel() {
    for (Command cmd : CoreCommand.INFO.commands()) {
      assertTrue(cmd.getCommandDescription() != null || cmd.isHidden());
      switch (cmd.getReleaseLevel()) {
        case ALPHA:
          assertTrue(cmd.getCommandName(), cmd.isHidden());
          break;
        case BETA:
        case GA:
          assertFalse(cmd.getCommandName(), cmd.isHidden());
          break;
        default:
          break;
      }
    }
  }

  public void testSpecificLicenceFlags() {
    assertEquals("enable_rtg", ToolsCommand.VERSION.getLicenceKeyName());
  }
}
