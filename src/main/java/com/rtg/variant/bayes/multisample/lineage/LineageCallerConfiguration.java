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
package com.rtg.variant.bayes.multisample.lineage;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.reference.Sex;
import com.rtg.reference.SexMemo;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.vcf.VariantStatistics;
import com.rtg.variant.bayes.complex.DenovoChecker;
import com.rtg.variant.bayes.complex.LineageDenovoChecker;
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

/**
 */
public final class LineageCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public LineageCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using Lineage caller");

      final String[] outputSampleNames = SamUtils.getSampleNames(params.uberHeader());
      final List<String> genomes = Arrays.asList(outputSampleNames);
      final int[] parents = new int[outputSampleNames.length];

      final GenomeRelationships genomeRelationships = params.genomeRelationships();
      if (outputSampleNames.length != genomeRelationships.genomes().length) {
        throw new NoTalkbackSlimException("Number of samples in pedigree and number of output genomes don't match");
      }
      final Lineage.LineageBuilder builder = new Lineage.LineageBuilder();
      final boolean noRelationships = genomeRelationships.relationships(Relationship.RelationshipType.PARENT_CHILD).length < 1
        && genomeRelationships.relationships(Relationship.RelationshipType.ORIGINAL_DERIVED).length < 1;
      if (noRelationships) {
        throw new NoTalkbackSlimException("Lineage requires at least one original derived or parent child relationship");
      }
      for (Relationship r : genomeRelationships.relationships(Relationship.RelationshipType.ORIGINAL_DERIVED)) {
        final int childPos = genomes.indexOf(r.second());
        parents[childPos]++;
        final int parentPos = genomes.indexOf(r.first());
        builder.add(parentPos, childPos);
      }
      for (Relationship r : genomeRelationships.relationships(Relationship.RelationshipType.PARENT_CHILD)) {
        final int childPos = genomes.indexOf(r.second());
        parents[childPos]++;
        final int parentPos = genomes.indexOf(r.first());
        builder.add(parentPos, childPos);
      }
      for (int i = 0; i < parents.length; ++i) {
        if (parents[i] > 1) {
          throw new NoTalkbackSlimException("Lineage requires at most one parent per individual, sample '" + outputSampleNames[i] + "' had " + parents[i]);
        }
      }

      if (outputSampleNames.length == 0) {
        throw new NoTalkbackSlimException("VCF output for lineage calling needs SAM headers with sample names");
      }

      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
      final PopulationHwHypothesesCreator ssp;
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      builder.deNovoPriorDefault(params.genomePriors().denovoRef()).params(params).coverage(false);
      final Lineage caller = builder.create();
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>(genomes.size());
      final SexMemo sexMemo = Utils.createSexMemo(params);
      Sex lineageSex = null;
      for (String genome : genomes) {
        final Sex sex = genomeRelationships.getSex(genome);
        if (lineageSex == null) {
          lineageSex = sex;
        } else if (lineageSex != sex) {
          throw new NoTalkbackSlimException("Lineage requires all samples to have the same sex");
        }
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, sex, sexMemo));
      }
      return new LineageCallerConfiguration(caller, outputSampleNames, individualFactories, chooser, haploid, diploid, ssp);
    }
  }

  private final LineageDenovoChecker mChecker;

  private LineageCallerConfiguration(Lineage jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
    mChecker = new LineageDenovoChecker(jointCaller.lineageLookup());
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP, VcfFormatField.SCONT));
    return f;
  }

  @Override
  public DenovoChecker getDenovoChecker() {
    return mChecker;
  }

  @Override
  public boolean handlesPloidy(final Ploidy ploidy) {
    return ploidy != Ploidy.POLYPLOID;
  }
}
