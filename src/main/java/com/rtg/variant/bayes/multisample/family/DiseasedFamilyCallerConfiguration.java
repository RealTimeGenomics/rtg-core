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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import com.rtg.reference.Sex;
import com.rtg.reference.SexMemo;
import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.PedigreeException;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;
import com.rtg.variant.bayes.multisample.IndividualSampleFactory;
import com.rtg.variant.bayes.multisample.JointCallerConfigurator;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.MultisampleUtils;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;
import com.rtg.vcf.VariantStatistics;

/**
 */
public final class DiseasedFamilyCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public DiseasedFamilyCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using disease caller");

      final GenomeRelationships genomeRelationships = params.genomeRelationships();
      final Family family;
      try {
        family = Family.getFamily(genomeRelationships);
      } catch (PedigreeException e) {
        throw new NoTalkbackSlimException(e.getMessage());
      }
      final String fatherName = family.getFather();
      final String motherName = family.getMother();
      final String[] genomes = new String[family.getChildren().length + Family.FIRST_CHILD_INDEX];
      final Sex[] sexes = new Sex[genomes.length];
      final List<String> childNames = new ArrayList<>(Arrays.asList(family.getChildren()));
      final String[] outputGenomes = SamUtils.getSampleNames(params.uberHeader());
      if (outputGenomes.length != genomes.length) {
        throw new NoTalkbackSlimException("Exactly " + genomes.length + " sample names expected in mappings");
      }
      genomes[Family.FATHER_INDEX] = fatherName;
      genomes[Family.MOTHER_INDEX] = motherName;
      sexes[Family.FATHER_INDEX] = genomeRelationships.getSex(fatherName);
      sexes[Family.MOTHER_INDEX] = genomeRelationships.getSex(motherName);
      Diagnostic.developerLog("First parent: " + fatherName + " sex: " + sexes[Family.FATHER_INDEX]);
      Diagnostic.developerLog("Second parent: " + motherName + " sex: " + sexes[Family.MOTHER_INDEX]);
      for (int i = 0; i < childNames.size(); ++i) {
        final String childName = childNames.get(i);
        genomes[Family.FIRST_CHILD_INDEX + i] = childName;
        sexes[Family.FIRST_CHILD_INDEX + i] = genomeRelationships.getSex(childName);
        Diagnostic.developerLog("Child " + i + ": " + childName + " sex: " + sexes[Family.FIRST_CHILD_INDEX + i]);
      }
      final List<String> familyNames = Arrays.asList(genomes);
      for (String mapName : outputGenomes) {
        if (!familyNames.contains(mapName)) {
          throw new NoTalkbackSlimException("Unexpected sample name in mappings: " + mapName);
        }
      }
      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
      final DiseasedFamilyCaller jointCaller = new DiseasedFamilyCaller(params, family, params.noDiseasePrior());
      final PopulationHwHypothesesCreator ssp;
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      final SexMemo sexMemo = Utils.createSexMemo(params);
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>(genomes.length);
      for (int i = 0; i < genomes.length; ++i) {
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, sexes[i], sexMemo));
      }
      return new DiseasedFamilyCallerConfiguration(jointCaller, genomes, individualFactories, chooser, haploid, diploid, ssp);
    }
  }

  private DiseasedFamilyCallerConfiguration(MultisampleJointCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP));
    f.addExtraInfoFields(EnumSet.of(VcfInfoField.DISEASE, VcfInfoField.RDS));
    return f;
  }
}
