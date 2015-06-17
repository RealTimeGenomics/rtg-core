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
package com.rtg.variant.bayes.multisample.family;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Set;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.relation.GenomeRelationships;
import com.rtg.usage.UsageMetric;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;

/**
 * Entry point to mendelian SNP calling.
 */
public class FamilyCli extends AbstractMultisampleCli {

  // Either this ...
  private static final String PEDIGREE_FLAG = "pedigree";

  // or these
  private static final String FATHER_FLAG = "father";
  private static final String MOTHER_FLAG = "mother";
  private static final String SON_FLAG = "son";
  private static final String DAUGHTER_FLAG = "daughter";

  private static final String MISMATCHED_PARAMS_ERROR1 = "Cannot use --" + PEDIGREE_FLAG + " in conjunction with --" + FATHER_FLAG + ", --" + MOTHER_FLAG + ", --" + SON_FLAG + ", and --" + DAUGHTER_FLAG + " flags";
  private static final String MISMATCHED_PARAMS_ERROR2 = "Must set --" + FATHER_FLAG + " and --" + MOTHER_FLAG + " flags and at least one of --" + SON_FLAG + " or --" + DAUGHTER_FLAG + " flags";

  static class FamilyValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (flags.isSet(PEDIGREE_FLAG)) {
        if (flags.isSet(FATHER_FLAG) || flags.isSet(MOTHER_FLAG) || flags.isSet(SON_FLAG) || flags.isSet(DAUGHTER_FLAG)) {
          flags.error(MISMATCHED_PARAMS_ERROR1);
          return false;
        }
      } else if (!flags.isSet(FATHER_FLAG) || !flags.isSet(MOTHER_FLAG) || !(flags.isSet(SON_FLAG) || flags.isSet(DAUGHTER_FLAG))) {
        flags.error(MISMATCHED_PARAMS_ERROR2);
        return false;
      } else {
        final String father = (String) flags.getValue(FATHER_FLAG);
        final String mother = (String) flags.getValue(MOTHER_FLAG);
        if (father.equals(mother)) {
          flags.error("Father and mother must be different samples");
          return false;
        } else {
          final Set<String> sons = new HashSet<>();
          final Set<String> daughters = new HashSet<>();
          if (flags.isSet(SON_FLAG)) {
            for (final Object son : flags.getValues(SON_FLAG)) {
              if (!sons.add((String) son)) {
                flags.error("Individual sons must be different sample to other sons");
                return false;
              }
            }
          }
          if (flags.isSet(DAUGHTER_FLAG)) {
            for (final Object daughter : flags.getValues(DAUGHTER_FLAG)) {
              if (!daughters.add((String) daughter)) {
                flags.error("Individual daughters must be different sample to other daughters");
                return false;
              }
            }
          }
          if (sons.contains(father) || sons.contains(mother)) {
            flags.error("Son must be different sample to mother and father");
            return false;
          } else if (daughters.contains(father) || daughters.contains(mother)) {
            flags.error("Daughter must be different sample to mother and father");
            return false;
          } else {
            final Set<String> combined = new HashSet<>();
            combined.addAll(sons);
            combined.addAll(daughters);
            if (combined.size() != sons.size() + daughters.size()) {
              flags.error("Son and daughter samples must be different");
              return false;
            }
          }
        }
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "family";
  }

  @Override
  public String description() {
    return "call variants for a family following Mendelian inheritance";
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    CommonFlags.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    flags.setDescription("Performs a combined mendelian family variant analysis.");
    flags.setValidator(new FamilyValidator());
    final Flag inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);

    final Flag pedFlag = flags.registerOptional('p', PEDIGREE_FLAG, File.class, "file", "genome relationships PED file, if not specifying family via --" + FATHER_FLAG + ", --" + MOTHER_FLAG + ", etc").setCategory(INPUT_OUTPUT);
    final Flag fatherFlag = flags.registerOptional(FATHER_FLAG, String.class, "string", "sample identifier used in read groups for father sample").setCategory(INPUT_OUTPUT);
    final Flag motherFlag = flags.registerOptional(MOTHER_FLAG, String.class, "string", "sample identifier used in read groups for for mother sample").setCategory(INPUT_OUTPUT);
    flags.registerOptional(SON_FLAG, String.class, "string", "sample identifier used in read groups for a son sample").setMaxCount(Integer.MAX_VALUE).setCategory(INPUT_OUTPUT);
    flags.registerOptional(DAUGHTER_FLAG, String.class, "string", "sample identifier used in read groups for a daughter sample").setMaxCount(Integer.MAX_VALUE).setCategory(INPUT_OUTPUT);
    AbstractMultisampleCli.registerComplexPruningFlags(flags);
    flags.addRequiredSet(fatherFlag, motherFlag, inFlag);
    flags.addRequiredSet(fatherFlag, motherFlag, listFlag);
    flags.addRequiredSet(pedFlag, inFlag);
    flags.addRequiredSet(pedFlag, listFlag);
  }

  @Override
  protected GenomeRelationships grf() throws IOException {
    if (mFlags.isSet(PEDIGREE_FLAG)) {
      final File relfile = (File) mFlags.getValue(PEDIGREE_FLAG);
      return GenomeRelationships.loadGenomeRelationships(relfile);
    } else {
      assert mFlags.isSet(FATHER_FLAG);
      assert mFlags.isSet(MOTHER_FLAG);
      assert mFlags.isSet(SON_FLAG) || mFlags.isSet(DAUGHTER_FLAG);
      final GenomeRelationships ped = new GenomeRelationships();

      final String father = (String) mFlags.getValue(FATHER_FLAG);
      final String mother = (String) mFlags.getValue(MOTHER_FLAG);
      ped.addGenome(mother, GenomeRelationships.SEX_FEMALE);
      ped.addGenome(father, GenomeRelationships.SEX_MALE);

      if (mFlags.isSet(SON_FLAG)) {
        for (final Object oson : mFlags.getValues(SON_FLAG)) {
          final String son = (String) oson;
          ped.addGenome(son, GenomeRelationships.SEX_MALE);
          ped.addParentChild(father, son);
          ped.addParentChild(mother, son);
        }
      }
      if (mFlags.isSet(DAUGHTER_FLAG)) {
        for (final Object odaughter : mFlags.getValues(DAUGHTER_FLAG)) {
          final String daughter = (String) odaughter;
          ped.addGenome(daughter, GenomeRelationships.SEX_FEMALE);
          ped.addParentChild(father, daughter);
          ped.addParentChild(mother, daughter);
        }
      }
      return ped;
    }
  }

  @Override
  public ParamsTask<?, ?> task(final VariantParams params, final OutputStream out) throws IOException {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MultisampleTask<>(params, new FamilyCallerConfiguration.Configurator(), out,  getStatistics(params), usageMetric);
  }

}
