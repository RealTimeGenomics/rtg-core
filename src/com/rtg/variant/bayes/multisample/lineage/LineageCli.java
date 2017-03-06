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
package com.rtg.variant.bayes.multisample.lineage;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.relation.GenomeRelationships;
import com.rtg.usage.UsageMetric;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.VariantParams;
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;

/**
 */
public class LineageCli extends AbstractMultisampleCli {

  private static final String PEDIGREE_FLAG = "pedigree";

  @Override
  protected GenomeRelationships grf() throws IOException {
    final File relfile = (File) mFlags.getValue(PEDIGREE_FLAG);
    return GenomeRelationships.loadGenomeRelationships(relfile);
  }

  @Override
  public String moduleName() {
    return "lineage";
  }

  @Override
  public String description() {
    return "call de novo variants in a cell lineage";
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }
  static class LineageValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      final File pedigreeFile = (File) flags.getValue(PEDIGREE_FLAG);
      if (!pedigreeFile.exists() || pedigreeFile.isDirectory()) {
        flags.setParseMessage("Pedigree file isn't readable");
        return false;
      }
      return true;
    }
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    AvrUtils.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    flags.setDescription("Performs a combined cell lineage variant analysis.");
    flags.setValidator(new LineageValidator());
    final Flag<File> inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerRequired('p', PEDIGREE_FLAG, File.class, "file", "genome relationships PED file").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
    registerComplexPruningFlags(flags, true);
  }
  @Override
  public ParamsTask<?, ?> task(final VariantParams params, final OutputStream out) throws IOException {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MultisampleTask<>(params, new LineageCallerConfiguration.Configurator(), out,  getStatistics(params), usageMetric);
  }

  /**
   * Command line
   * @param args command line args
   */
  public static void main(String[] args) {
    new LineageCli().mainInit(args, System.out, System.err);
  }
}
