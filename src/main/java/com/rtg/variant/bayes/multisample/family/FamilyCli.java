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
package com.rtg.variant.bayes.multisample.family;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_LIST_FLAG;
import static com.rtg.launcher.CommonFlags.PEDIGREE_FLAG;
import static com.rtg.launcher.CommonFlags.STRING;
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
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;

/**
 * Entry point to mendelian SNP calling.
 */
public class FamilyCli extends AbstractMultisampleCli {

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
          flags.setParseMessage(MISMATCHED_PARAMS_ERROR1);
          return false;
        }
      } else if (!flags.isSet(FATHER_FLAG) || !flags.isSet(MOTHER_FLAG) || !(flags.isSet(SON_FLAG) || flags.isSet(DAUGHTER_FLAG))) {
        flags.setParseMessage(MISMATCHED_PARAMS_ERROR2);
        return false;
      } else {
        final String father = (String) flags.getValue(FATHER_FLAG);
        final String mother = (String) flags.getValue(MOTHER_FLAG);
        if (father.equals(mother)) {
          flags.setParseMessage("Father and mother must be different samples");
          return false;
        } else {
          final Set<String> sons = new HashSet<>();
          final Set<String> daughters = new HashSet<>();
          if (flags.isSet(SON_FLAG)) {
            for (final Object son : flags.getValues(SON_FLAG)) {
              if (!sons.add((String) son)) {
                flags.setParseMessage("Individual sons must be different sample to other sons");
                return false;
              }
            }
          }
          if (flags.isSet(DAUGHTER_FLAG)) {
            for (final Object daughter : flags.getValues(DAUGHTER_FLAG)) {
              if (!daughters.add((String) daughter)) {
                flags.setParseMessage("Individual daughters must be different sample to other daughters");
                return false;
              }
            }
          }
          if (sons.contains(father) || sons.contains(mother)) {
            flags.setParseMessage("Son must be different sample to mother and father");
            return false;
          } else if (daughters.contains(father) || daughters.contains(mother)) {
            flags.setParseMessage("Daughter must be different sample to mother and father");
            return false;
          } else {
            final Set<String> combined = new HashSet<>();
            combined.addAll(sons);
            combined.addAll(daughters);
            if (combined.size() != sons.size() + daughters.size()) {
              flags.setParseMessage("Son and daughter samples must be different");
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
    AvrUtils.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    flags.setDescription("Performs a combined mendelian family variant analysis.");
    flags.setValidator(new FamilyValidator());
    final Flag<File> inFlag = flags.registerRequired(File.class, FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', INPUT_LIST_FLAG, File.class, FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);

    final Flag<File> pedFlag = flags.registerOptional('p', PEDIGREE_FLAG, File.class, FILE, "genome relationships PED file, if not specifying family via --" + FATHER_FLAG + ", --" + MOTHER_FLAG + ", etc").setCategory(INPUT_OUTPUT);
    final Flag<String> fatherFlag = flags.registerOptional(FATHER_FLAG, String.class, STRING, "sample identifier used in read groups for father sample").setCategory(INPUT_OUTPUT);
    final Flag<String> motherFlag = flags.registerOptional(MOTHER_FLAG, String.class, STRING, "sample identifier used in read groups for for mother sample").setCategory(INPUT_OUTPUT);
    flags.registerOptional(SON_FLAG, String.class, STRING, "sample identifier used in read groups for a son sample").setMaxCount(Integer.MAX_VALUE).setCategory(INPUT_OUTPUT);
    flags.registerOptional(DAUGHTER_FLAG, String.class, STRING, "sample identifier used in read groups for a daughter sample").setMaxCount(Integer.MAX_VALUE).setCategory(INPUT_OUTPUT);
    registerComplexPruningFlags(flags, true);
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
