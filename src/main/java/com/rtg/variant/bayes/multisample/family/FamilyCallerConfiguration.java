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

import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;
import com.rtg.reference.SexMemo;
import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.PedigreeException;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.complex.DenovoChecker;
import com.rtg.variant.bayes.complex.MendelianDenovoChecker;
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
import com.rtg.vcf.ChildPhasingVcfAnnotator;
import com.rtg.vcf.VariantStatistics;

/**
 */
public final class FamilyCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public FamilyCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using Mendelian caller");

      final String[] outputGenomes = SamUtils.getSampleNames(params.uberHeader());
      if (outputGenomes.length == 0) {
        throw new NoTalkbackSlimException("VCF output for family calling needs SAM headers with sample names");
      }

      final List<String> tmp = Arrays.asList(outputGenomes);

      final GenomeRelationships genomeRelationships = params.genomeRelationships();
      final List<String> calledGenomes = new ArrayList<>();

      final Family family;
      try {
        family = Family.getFamily(genomeRelationships);
        for (String child : family.getChildren()) {
          if (genomeRelationships.getSex(child) == Sex.EITHER) {
            throw new NoTalkbackSlimException("The sex of child " + child + " was not specified");
          }
        }
      } catch (PedigreeException e) {
        throw new NoTalkbackSlimException("There was a problem with the family pedigree: " + e.getMessage());
      }
      Diagnostic.developerLog(family.toString());

      // First pass to prefill calledGenomes with only the family members we want the calls to be output for.
      for (String member : family.getMembers()) {
        if (tmp.contains(member)) {
          calledGenomes.add(member);
        }
      }

      final String[] newOutputGenomes = calledGenomes.toArray(new String[0]);
      if (newOutputGenomes.length != outputGenomes.length) {
        throw new NoTalkbackSlimException("Not all sample names are members of the family");
      }

      // Check family has enough called samples to be usable
      if (!Utils.isCallableAsFamily(calledGenomes, family)) {
        throw new NoTalkbackSlimException("Not enough family members have mapping data provided");
      }

      // Second pass will add any family members for which calls will be made but not output
      for (String member : family.getMembers()) {
        if (!calledGenomes.contains(member)) {
          calledGenomes.add(member);
        }
      }

      // Set sample ids
      family.setSampleIds(calledGenomes);

      if (!ReferenceGenome.hasReferenceFile(params.genome().reader())) {
        Diagnostic.warning("The reference SDF does not contain reference.txt configuration, assuming autosomal inheritance (calling on sex chromosomes will not be correct). For more information, see the user manual.");
      }

      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
      final PopulationHwHypothesesCreator ssp;
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      final FamilyCaller familyCaller = new FamilyCaller(params, family);
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>();
      final SexMemo sexMemo = Utils.createSexMemo(params);
      for (String genome : calledGenomes) {
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, genomeRelationships.getSex(genome), sexMemo));
      }
      final FamilyCallerConfiguration fc = new FamilyCallerConfiguration(familyCaller, newOutputGenomes, individualFactories, chooser, haploid, diploid, ssp, new ChildFamilyLookup(calledGenomes.size(), family));
      fc.getVcfAnnotators().add(new ChildPhasingVcfAnnotator(family));
      if (params.infoAnnotations().contains(VcfInfoField.SGP)) {
        fc.getVcfAnnotators().add(new SegregationVcfAnnotator(family));
      }
      return fc;
    }
  }

  private final DenovoChecker mDenovoCorrector;

  private FamilyCallerConfiguration(MultisampleJointCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp, ChildFamilyLookup family) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
    mDenovoCorrector = new MendelianDenovoChecker(family);
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP, VcfFormatField.SCONT));
    return f;
  }

  @Override
  public DenovoChecker getDenovoChecker() {
    return mDenovoCorrector;
  }

}
