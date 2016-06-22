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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.IOException;

import com.rtg.launcher.CommonFlags;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class TumorOnlyCli extends SomaticCli {
  private static final String ORIGINAL_NAME = "normal";
  protected static final String SAMPLE_FLAG = "sample";

  @Override
  protected GenomeRelationships grf() throws IOException {
      assert mFlags.isSet(SAMPLE_FLAG);

      final GenomeRelationships ped = new GenomeRelationships();
      final String original = ORIGINAL_NAME;
      final String derived = (String) mFlags.getValue(SAMPLE_FLAG);
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
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, TEMPLATE_FLAG)) {
        return false;
      }
      if (!flags.checkInRange(SOMATIC_FLAG, 0.0, false, 1.0, false)
        || !flags.checkInRange(LOH_FLAG, 0.0, true, 1.0, true)
        || !flags.checkInRange(CONTAMINATION_FLAG, 0.0, true, 1.0, false)
        || !flags.checkInRange(REVERSE_CONTAMINATION_FLAG, 0.0, true, 1.0, false)) {
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
    flags.setDescription("Performs a somatic variant analysis on a mixed tumor sample.");
    flags.registerRequired(SAMPLE_FLAG, String.class, "string", "sample identifier used in read groups for tumor sample").setCategory(INPUT_OUTPUT);
    final Flag contamination = flags.getFlag(CONTAMINATION_FLAG);
    contamination.setParameterDefault(0.75);
    requiredSet(flags);
    flags.setValidator(new MultigenomeValidator());
  }

  private void requiredSet(CFlags flags) {
    flags.addRequiredSet(flags.getAnonymousFlag(0));
    flags.addRequiredSet(flags.getFlag(CommonFlags.INPUT_LIST_FLAG));
  }

  @Override
  protected SomaticParamsBuilder getSomaticParamsBuilder() {
    return super.getSomaticParamsBuilder().includeGermlineVariants(true);
  }
}
