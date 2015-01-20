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
import java.util.Collections;
import java.util.HashSet;

import com.rtg.launcher.GlobalFlags;
import com.rtg.util.TestUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests corresponding class
 */
public class CoreCommandTest extends TestCase {

  @Override
  public void setUp() {
    GlobalFlags.resetAccessedStatus();
  }

  public void testEnum() {
    // Enum, not display order
    TestUtils.testEnum(CoreCommand.class, "[FORMAT, CG2SDF, SDF2FASTA, SDF2FASTQ, SDF2QUALA, MAP, MAPF, CGMAP, ZBUILD, ZMAP, ZMAPF, MAPX, ASSEMBLE, ADDPACBIO, "
      + "SVPREP, SV, DISCORD, COVERAGE, SNP, MULTISNP, FAMILY, SOMATIC, POPULATION, LINEAGE, AVRBUILD, AVRPREDICT, CNV, CALIBRATE, SPECIES, METASNP, SIMILARITY, COMPOSITION-META-PIPELINE, FUNCTIONAL-META-PIPELINE, COMPOSITION-FUNCTIONAL-META-PIPELINE, GENOMESIM, CGSIM, READSIM, READSIMEVAL, SNPSIM, POPSIM, SAMPLESIM, CHILDSIM, DENOVOSIM, SAMPLEREPLAY, VCFEVAL, CNVSIM, CNVSIMEVAL, BGZIP, INDEX, EXTRACT, AVIEW, SDFSTATS, SDFSPLIT, SDFSUBSET, SDFSUBSEQ, SAM2BAM, SAMMERGE, SAMSTATS, SAMRENAME, MAPXRENAME, CHRSTATS, VCFFILTER, VCFANNOTATE, SNPINTERSECT, NCBI2TAX, TAXFILTER, TAXSTATS, MENDELIAN, PHASINGSEARCH, PHASINGEVAL, VCFSTATS, VCFMERGE, VCFSUBSET, PEDFILTER, PEDSTATS, AVRSTATS, ROCPLOT, USAGESERVER, VERSION, LICENSE, HELP]");
    final HashSet<Command> displayCommands = new HashSet<>();
    Collections.addAll(displayCommands, CoreCommand.INFO.commands());
    for (CoreCommand mod : CoreCommand.values()) {
      assertTrue(displayCommands.contains(mod.module()));
    }
  }

  public void testUsage() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    CoreCommand.FORMAT.module().mainInit(new String[]{"-h"}, baos, null);
    baos.flush();
    TestUtils.containsAll(baos.toString().replaceAll("\\s+", " "), "Usage: rtg format [OPTION]... -o SDF FILE+",
        "-h, --help print help on command-line flag usage");
    assertEquals("FORMAT", CoreCommand.FORMAT.module().getCommandName());
    assertEquals("FORMAT", CoreCommand.FORMAT.toString());
    assertEquals(CommandCategory.FORMAT, CoreCommand.FORMAT.module().getCategory());
    assertEquals(CoreCommand.FORMAT, CoreCommand.valueOf("FORMAT"));
  }

  public void testModules() {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(baos);
    assertEquals(0, CoreCommand.VERSION.module().mainInit(new String[0], baos, ps));
    assertEquals(0, CoreCommand.LICENSE.module().mainInit(new String[0], baos, ps));
    assertEquals(0, CoreCommand.HELP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CG2SDF.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CGMAP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CGSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNV.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CALIBRATE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MENDELIAN.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.COVERAGE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CHRSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.FORMAT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.GENOMESIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPX.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.COMPOSITIONMETA.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.FUNCTIONALMETA.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.METAGENOMICS.module().mainInit(new String[0], baos, ps));

    assertEquals(1, CoreCommand.SIMILARITY.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPF.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.READSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAM2BAM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMRENAME.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMMERGE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MAPXRENAME.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDF2FASTA.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDF2FASTQ.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDF2QUALA.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDFSUBSET.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDFSUBSEQ.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDFSPLIT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SDFSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFFILTER.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFANNOTATE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SNPINTERSECT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.TAXFILTER.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.NCBI2TAX.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.TAXSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFMERGE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SNPSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.POPSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMPLESIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CHILDSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.DENOVOSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SAMPLEREPLAY.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFEVAL.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SINGLETON.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.READSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSIM.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.READSIMEVAL.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.CNVSIMEVAL.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.ASSEMBLE.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.ADDPACBIO.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.INDEX.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVIEW.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SPECIES.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.METASNP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SVPREP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.SV.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.DISCORD.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_VARIANT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_FAMILY.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_SOMATIC.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.MULTI_POPULATION.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.BGZIP.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.EXTRACT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRBUILD.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRPREDICT.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.AVRSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.PEDFILTER.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.PEDSTATS.module().mainInit(new String[0], baos, ps));
    assertEquals(1, CoreCommand.VCFSUBSET.module().mainInit(new String[0], baos, ps));
  }

  public void testModuleHelp() {
    final String[] helpArgs = {"--help"};
    for (final CoreCommand module : CoreCommand.values()) {
      if ((module == CoreCommand.VERSION) || (module == CoreCommand.LICENSE) || (module == CoreCommand.HELP)) {
        continue; // These modules have no help
      }
      GlobalFlags.resetAccessedStatus();
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      final PrintStream ps = new PrintStream(baos);
      module.module().mainInit(helpArgs, baos, ps);
      assertTrue(module.module().getCommandName() + " produces no help output: " + baos.toString(), baos.toString().contains("Usage"));
    }
  }

  public void testReleaseLevel() {
    for (CoreCommand cmd : CoreCommand.values()) {
      final Command mod = cmd.module();
      switch (mod.getReleaseLevel()) {
        case ALPHA:
          assertTrue(cmd.name(), mod.isHidden());
          break;
        case BETA:
        case GA:
          assertFalse(cmd.name(), mod.isHidden());
          break;
        default:
          break;
      }
    }
  }

  public void testSpecificLicenceFlags() {
    assertEquals("enable_rtg", CoreCommand.VERSION.module().getLicenceKeyName());
  }

  public static Test suite() {
    return new TestSuite(CoreCommandTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(CoreCommandTest.class));
  }
}
