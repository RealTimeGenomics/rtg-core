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

package com.rtg.variant.bayes.multisample.statistics;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamFilterOptions;
import com.rtg.usage.UsageMetric;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;
import com.rtg.vcf.VariantStatistics;

/**
 */
public class VstatsCli extends AbstractMultisampleCli {

  private static final String SEX_FLAG = "sex";
  private static final String PLOIDY_FLAG = "ploidy";

  @Override
  protected GenomeRelationships grf() {
    return null;
  }

  @Override
  public String moduleName() {
    return "vstats";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  @TestClass("com.rtg.variant.bayes.multisample.singleton.SingletonValidatorTest")
  static class SingletonValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, CommonFlags.TEMPLATE_FLAG)) {
        return false;
      }
      return true;
    }
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    flags.setDescription("Calls sequence variants, such as single nucleotide polymorphisms (SNPs), multi-nucleotide polymorphisms (MNPs) and Indels, from a set of alignments reported in genome-sorted SAM/BAM position.");
    flags.setValidator(new SingletonValidator());
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setMinCount(0);
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(PLOIDY_FLAG, ReferencePloidy.class, CommonFlags.STRING, "ploidy to use", ReferencePloidy.AUTO).setCategory(SENSITIVITY_TUNING);

    SamFilterOptions.registerMaxHitsFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASMatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASUnmatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerExcludeMatedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);
    flags.registerOptional(SEX_FLAG, Sex.class, "sex", "sex of individual", Sex.EITHER).setCategory(SENSITIVITY_TUNING);

    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  @Override
  protected VariantParamsBuilder makeParamsBuilder() throws IOException {
    final VariantParamsBuilder builder = super.makeParamsBuilder();
    builder.sex((Sex) mFlags.getValue(SEX_FLAG));
    builder.ploidy((ReferencePloidy) mFlags.getValue(PLOIDY_FLAG));
    builder.noComplexCalls(true);
    return builder;
  }

  @Override
  public VariantParams makeParams() throws IOException {
    return super.makeParams();
  }

  @Override
  protected ParamsTask<?, ?> task(final VariantParams params, final OutputStream out) throws IOException {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MultisampleTask<>(params, new VstatsCallerConfiguration.Configurator(), out, new VariantStatistics(params.directory()), usageMetric);
  }

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new VstatsCli().mainExit(args);
  }
}
