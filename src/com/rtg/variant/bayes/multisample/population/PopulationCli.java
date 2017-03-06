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

package com.rtg.variant.bayes.multisample.population;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.usage.UsageMetric;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.GenomeConnectivity;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;

/**
 */
public class PopulationCli extends AbstractMultisampleCli {

  private static final String PEDIGREE_FLAG = "pedigree";

  private static final String IMPUTE_FLAG = "impute";

  private static final String REMOVE_RELATIONSHIPS_FLAG = "Xremove-relationships";

  private static final String PROPAGATING_FLAG = "Xprop-priors";

  private static final String MAX_EM_ITERATIONS_FLAG = "Xmax-em-iterations";

  private static final String PEDIGREE_CONNECTIVITY_FLAG = "pedigree-connectivity";

  static class MultisampleValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (flags.isSet(MAX_EM_ITERATIONS_FLAG) && ((Integer) flags.getValue(MAX_EM_ITERATIONS_FLAG) < 0)) {
        flags.error("Invalid value for maximum number of EM iterations, minimum is 0");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "population";
  }

  @Override
  public String description() {
    return "call variants for multiple potentially-related individuals";
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    AvrUtils.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    flags.setDescription("Performs a multiple sample variant analysis of many, potentially related, genomes.");
    flags.setValidator(new MultisampleValidator());
    final Flag<File> inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads")
        .setCategory(INPUT_OUTPUT)
        .setMinCount(0)
        .setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerRequired('p', PEDIGREE_FLAG, File.class, "file", "genome relationships PED file").setCategory(INPUT_OUTPUT);
    flags.registerOptional(IMPUTE_FLAG, String.class, "STRING", "name of sample absent from mappings to impute genotype for")
    .setCategory(REPORTING)
    .setMaxCount(Integer.MAX_VALUE).enableCsv();
    registerComplexPruningFlags(flags, true);
    flags.registerOptional(REMOVE_RELATIONSHIPS_FLAG, "if set, remove all relationship information from pedigree.").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(PROPAGATING_FLAG, "true to use propagating priors algorithm for pedigrees").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(MAX_EM_ITERATIONS_FLAG, Integer.class, "INT", "maximum number of EM iterations. 0 to disable EM", EmAlgorithm.DEFAULT_MAX_ITERATIONS).setCategory(REPORTING);

    final ArrayList<String> connectivties = new ArrayList<>();
    connectivties.add("auto");
    for (GenomeConnectivity conn : GenomeConnectivity.values()) {
      connectivties.add(conn.toString().toLowerCase(Locale.getDefault()));
    }
    flags.registerOptional(PEDIGREE_CONNECTIVITY_FLAG, String.class, "STRING", "sets mode of operation based on how well connected the pedigree is", "auto").setCategory(SENSITIVITY_TUNING).setParameterRange(connectivties);

    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  @Override
  protected GenomeRelationships grf() throws IOException {
    final File relfile = (File) mFlags.getValue(PEDIGREE_FLAG);
    GenomeRelationships gr = GenomeRelationships.loadGenomeRelationships(relfile);
    if (mFlags.isSet(REMOVE_RELATIONSHIPS_FLAG)) {
      gr = gr.filterByRelationships(new Relationship.NotFilter(new Relationship.RelationshipTypeFilter(Relationship.RelationshipType.PARENT_CHILD)));
    }
    return gr;
  }

  @Override
  protected VariantParamsBuilder makeParamsBuilder() throws InvalidParamsException, IOException {
    final VariantParamsBuilder builder = super.makeParamsBuilder();
    if (mFlags.isSet(IMPUTE_FLAG)) {
      final List<String> impute = new ArrayList<>();
      for (final Object o : mFlags.getValues(IMPUTE_FLAG)) {
        impute.add((String) o);
      }
      builder.imputedSamples(impute);
    }
    final boolean propagating = mFlags.isSet(PROPAGATING_FLAG);
    if (propagating) {
      builder.usePropagatingPriors(propagating);
      builder.maxEmIterations(0);
    } else if (mFlags.isSet(MAX_EM_ITERATIONS_FLAG)) {
      builder.maxEmIterations((Integer) mFlags.getValue(MAX_EM_ITERATIONS_FLAG));
    }

    final String conn = ((String) mFlags.getValue(PEDIGREE_CONNECTIVITY_FLAG)).toUpperCase(Locale.getDefault());
    if (!"AUTO".equals(conn)) {
      builder.genomeConnectivity(GenomeConnectivity.valueOf(conn));
    }

    return builder;
  }

  @Override
  protected ParamsTask<?, ?> task(final VariantParams params, final OutputStream out) throws IOException {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MultisampleTask<>(params, new PopulationCallerConfiguration.Configurator(), out,  getStatistics(params), usageMetric);
  }

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new PopulationCli().mainInit(args, System.out, System.err);
  }
}
