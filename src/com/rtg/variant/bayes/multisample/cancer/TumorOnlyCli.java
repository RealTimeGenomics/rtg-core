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

package com.rtg.variant.bayes.multisample.cancer;

import java.io.IOException;

import com.rtg.launcher.CommonFlags;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class TumorOnlyCli extends SomaticCli {
  private static final String ORIGINAL_NAME = "normal";

  @Override
  protected GenomeRelationships grf() throws IOException {
      assert mFlags.isSet(DERIVED_FLAG);
      assert mFlags.isSet(CONTAMINATION_FLAG);

      final GenomeRelationships ped = new GenomeRelationships();
      final String original = ORIGINAL_NAME;
      final String derived = (String) mFlags.getValue(DERIVED_FLAG);
      final Sex sex = (Sex) mFlags.getValue(SEX_FLAG);
      ped.addGenome(original, sex);
      ped.addGenome(derived, sex);
      final Relationship relationship = ped.addRelationship(Relationship.RelationshipType.ORIGINAL_DERIVED, original, derived);
      relationship.setProperty(Relationship.CONTAMINATION, String.valueOf(mFlags.getValue(CONTAMINATION_FLAG)));
      relationship.setProperty(Relationship.REVERSE_CONTAMINATION, String.valueOf(mFlags.getValue(REVERSE_CONTAMINATION_FLAG)));
      return ped;
  }

  private static class MultigenomeValidator implements Validator {

    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (!flags.isSet(DERIVED_FLAG) || !flags.isSet(CONTAMINATION_FLAG)) {
        flags.error("both --" + DERIVED_FLAG + " and --" + CONTAMINATION_FLAG + " must be specified");
        return false;
      }
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, TEMPLATE_FLAG)) {
        return false;
      }
      final double som = (Double) flags.getValue(SOMATIC_FLAG);
      if (som <= 0 || som >= 1) {
        flags.error("--" + SOMATIC_FLAG + " should be a probability 0<s<1");
        return false;
      }
      final double loh = (Double) flags.getValue(LOH_FLAG);
      if (loh < 0 || loh > 1) {
        flags.error("--" + LOH_FLAG + " should be a probability 0<=p<=1");
        return false;
      }
      final double reverse = (Double) flags.getValue(REVERSE_CONTAMINATION_FLAG);
      if (reverse < 0 || reverse > 1) {
        flags.error("--" + REVERSE_CONTAMINATION_FLAG + " should be a probability 0<p<1");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "tumoronly";
  }

  @Override
  public String description() {
    return "call variants for a mixed tumor/normal sample";
  }

  @Override
  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    CommonFlags.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    commonSomaticFlags(flags);
    flags.setDescription("Performs a somatic variant analysis on a mixed tumor/normal sample.");
    requiredSet(flags);
    flags.setValidator(new MultigenomeValidator());
  }
  private void requiredSet(CFlags flags) {
    flags.addRequiredSet(flags.getFlag(DERIVED_FLAG), flags.getFlag(CONTAMINATION_FLAG), flags.getAnonymousFlag(0));
    flags.addRequiredSet(flags.getFlag(DERIVED_FLAG), flags.getFlag(CONTAMINATION_FLAG), flags.getFlag(CommonFlags.INPUT_LIST_FLAG));
  }

  @Override
  protected SomaticParamsBuilder getSomaticParamsBuilder() {
    return super.getSomaticParamsBuilder().includeGermlineVariants(true);
  }
}
