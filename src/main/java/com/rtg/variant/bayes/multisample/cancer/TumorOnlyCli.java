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

package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import com.rtg.launcher.CommonFlags;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 * Somatic calling entry point for calling with a single sample.
 */
public class TumorOnlyCli extends SomaticCli {
  private static final String ORIGINAL_NAME = "normal";
  protected static final String SAMPLE_FLAG = "sample";

  @Override
  protected GenomeRelationships grf() {
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
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, CommonFlags.TEMPLATE_FLAG)) {
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
  @SuppressWarnings("unchecked")
  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    AvrUtils.initAvrModel(flags, false, SOMATIC_MODEL_DEFAULT);
    CommonFlags.initMinAvrScore(flags);
    commonSomaticFlags(flags);
    flags.setDescription("Performs a somatic variant analysis on a mixed tumor sample where no separate normal sample is available.");
    flags.registerRequired(SAMPLE_FLAG, String.class, CommonFlags.STRING, "sample identifier used in read groups for tumor sample").setCategory(INPUT_OUTPUT);
    final Flag<Double> contamination = (Flag<Double>) flags.getFlag(CONTAMINATION_FLAG);
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
