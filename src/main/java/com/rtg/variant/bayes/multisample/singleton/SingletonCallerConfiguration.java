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

package com.rtg.variant.bayes.multisample.singleton;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import com.rtg.reference.Sex;
import com.rtg.reference.SexMemo;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;
import com.rtg.variant.bayes.multisample.IndividualSampleFactory;
import com.rtg.variant.bayes.multisample.JointCallerConfigurator;
import com.rtg.variant.bayes.multisample.MultisampleUtils;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VariantStatistics;

/**
 */
public final class SingletonCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public SingletonCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using singleton caller");

      final String[] outputSampleNames = SamUtils.getSampleNames(params.uberHeader());
      if (outputSampleNames.length > 1) {
        throw new NoTalkbackSlimException("Multiple samples detected in read group headers, but this command models only one genome");
      }

      final String sampleName = (outputSampleNames.length == 0) ? VariantOutputVcfFormatter.DEFAULT_SAMPLE : outputSampleNames[0];

      final GenomeRelationships genomeRelationships = params.genomeRelationships();
      final Sex sampleSex;
      if (genomeRelationships != null) {
        if (outputSampleNames.length == 0) {
          throw new NoTalkbackSlimException("Pedigree was supplied but mappings do not contain a sample name");
        }
        if (!Arrays.asList(genomeRelationships.genomes()).contains(sampleName)) {
          throw new NoTalkbackSlimException("Mapping sample " + sampleName + " is not contained in the pedigree");
        }
        sampleSex = genomeRelationships.getSex(sampleName);
      } else {
        sampleSex = params.sex();
      }

      final String[] genomeNames = {sampleName};
      final SingletonCaller singletonCaller = new SingletonCaller(params);
      final PopulationHwHypothesesCreator ssp;
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>();
      final SexMemo sexMemo = Utils.createSexMemo(params);
      individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, sampleSex, sexMemo));
      return new SingletonCallerConfiguration(singletonCaller, genomeNames, individualFactories, chooser, haploid, diploid, ssp);
    }
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    if (params.minVariantAllelicDepth() > 0 || params.minVariantAllelicFraction() > 0.0) {
      f.addExtraFormatFields(EnumSet.of(VcfFormatField.VAF));
    }
    return f;
  }

  SingletonCallerConfiguration(SingletonCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
  }
}
