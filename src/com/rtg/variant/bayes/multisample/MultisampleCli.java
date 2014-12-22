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
package com.rtg.variant.bayes.multisample;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.CommonFlags;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.VariantParamsBuilder;

/**
 */
public class MultisampleCli extends AbstractMultisampleCli {

  private static final String MODULE_NAME = "multisnp";

  private static final String NO_DISEASE_FLAG = "Xno-disease-prior";
  private static final String RELATIONS_FLAG = "relations";


  static class MultisampleValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      final double dp = (Double) flags.getValue(NO_DISEASE_FLAG);
      if (dp <= 0 || dp >= 1) {
        flags.error("--" + NO_DISEASE_FLAG + " should be a probability 0<s<1");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    flags.setDescription("Entry point for alternative multiple sample variant analysis methods.");
    flags.setValidator(new MultisampleValidator());
    addSamFileFlags(flags);
    flags.registerRequired('r', RELATIONS_FLAG, File.class, "file", "relationships file").setCategory(INPUT_OUTPUT);
    flags.registerOptional('d', NO_DISEASE_FLAG, Double.class, "float", "prior probability of no disease per position", 0.95).setCategory(SENSITIVITY_TUNING);

    SamFilterOptions.registerMaxHitsFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASMatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASUnmatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerExcludeMatedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);

  }

  /**
   * Sam file flags (input/file list)
   * @param flags object to add flag to
   */
  public static void addSamFileFlags(CFlags flags) {
    final Flag inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  @Override
  protected GenomeRelationships grf() throws IOException {
    final File relfile = (File) mFlags.getValue(RELATIONS_FLAG);
    return GenomeRelationships.loadGenomeRelationships(relfile);
  }

  @Override
  protected VariantParamsBuilder makeParamsBuilder() throws InvalidParamsException, IOException {
    final VariantParamsBuilder builder = super.makeParamsBuilder();
    builder.noDiseasePrior((Double) mFlags.getValue(NO_DISEASE_FLAG));
    return builder;
  }
}
