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
package com.rtg.variant.bayes.multisample.population;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.SexMemo;
import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.MultiFamilyOrdering;
import com.rtg.relation.PedigreeException;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.GenomeConnectivity;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.complex.DenovoChecker;
import com.rtg.variant.bayes.complex.MendelianDenovoChecker;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;
import com.rtg.variant.bayes.multisample.IndividualSampleFactory;
import com.rtg.variant.bayes.multisample.JointCallerConfigurator;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.MultisampleUtils;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.family.FamilyCaller;
import com.rtg.variant.bayes.multisample.forwardbackward.FamilyCallerFB;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.ChildPhasingVcfAnnotator;
import com.rtg.vcf.VariantStatistics;
import com.rtg.vcf.VcfAnnotator;

/**
 */
public final class PopulationCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public PopulationCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using Population caller");
      final List<String> mappedGenomes = Arrays.asList(SamUtils.getSampleNames(params.uberHeader()));
      final GenomeRelationships genomeRelationships = params.genomeRelationships();
      final List<String> pedGenomes = Arrays.asList(genomeRelationships.genomes());
      final List<String> outputGenomes = new ArrayList<>(mappedGenomes.size() + params.imputedSamples().size());

      // Validate supplied names against what we're using from the relationships
      if (mappedGenomes.isEmpty()) {
        throw new NoTalkbackSlimException("VCF output for Multigenome needs SAM headers with sample names");
      }
      for (String genome : mappedGenomes) {
        if (!pedGenomes.contains(genome)) {
          throw new NoTalkbackSlimException("SAM sample " + genome + " is not contained in the pedigree");
        }
        outputGenomes.add(genome);
      }

      for (String genome : params.imputedSamples()) {
        if (!pedGenomes.contains(genome)) {
          throw new NoTalkbackSlimException("Imputed sample " + genome + " is not contained in the pedigree");
        }
        if (mappedGenomes.contains(genome)) {
          throw new NoTalkbackSlimException("Imputed sample " + genome + " is already contained in SAM samples");
        }
        outputGenomes.add(genome);
      }

      GenomeConnectivity gc = params.genomeConnectivity();
      if (gc == null) {
        gc = GenomeConnectivity.getConnectivity(outputGenomes, genomeRelationships);
        Diagnostic.userLog("Automatically determined genome connectivity to be: " + gc);
      }
      Diagnostic.developerLog("Genome connectivity: " + gc);
      // if dense connectivity update EM iterations to 0, unless they have been set elsewhere
      final VariantParams newParams;
      if (gc == GenomeConnectivity.DENSE && params.maxEmIterations() < 0) {
        newParams = params.cloneBuilder().maxEmIterations(0).create();
      } else {
        newParams = params;
      }

      final List<String> calledGenomes = new ArrayList<>(outputGenomes);
      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(newParams);
      final PopulationHwHypothesesCreator ssp;
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      MultisampleJointScorer familyCaller = null;
      Family[] famArray = null;
      final List<VcfAnnotator> annot = new ArrayList<>();
      try {
        // Only include families where at least two members are in the output genomes, at least one of which is a child
        final Set<Family> initialFamilies = Family.getFamilies(genomeRelationships, false, null);
        final Set<Family> families = new LinkedHashSet<>();
        for (Family family : initialFamilies) {
          if (Utils.isCallableAsFamily(outputGenomes, family)) {
            families.add(family);
          } else {
            Diagnostic.developerLog("Excluded family " + family + " due to not enough members included in output");
          }
        }

        if (families.size() > 0) {
          if (!ReferenceGenome.hasReferenceFile(params.genome().reader())) { // This isn't such a problem when no families are involved
            Diagnostic.warning("The reference SDF does not contain reference.txt configuration, assuming autosomal inheritance (calling on sex chromosomes will not be correct). For more information, see the user manual.");
          }

          // Ensure all family members are in called genomes
          for (Family family : families) {
            for (String member : family.getMembers()) {
              if (!calledGenomes.contains(member)) {
                calledGenomes.add(member);
              }
            }
          }
          // Set sample ID mappings
          for (Family family : families) {
            family.setSampleIds(calledGenomes);
          }
          Diagnostic.userLog("Identified " + families.size() + " usable families");
          Diagnostic.userLog("Families: " + StringUtils.LS + families);
          final List<Family> orderedFamilies = MultiFamilyOrdering.orderFamiliesAndSetMates(families);
          famArray = orderedFamilies.toArray(new Family[0]);
          familyCaller = newParams.usePropagatingPriors() ? new FamilyCallerFB(newParams, famArray) : new FamilyCaller(newParams, famArray);
          annot.add(new ChildPhasingVcfAnnotator(orderedFamilies));
        }
      } catch (PedigreeException e) {
        throw new NoTalkbackSlimException("There was a problem with the supplied pedigree: " + e.getMessage());
      }
      final PopulationCaller popCaller = new PopulationCaller(newParams, familyCaller);
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>(calledGenomes.size());
      final SexMemo sexMemo = Utils.createSexMemo(newParams);
      for (String genome : calledGenomes) {
        Diagnostic.userLog("Sample: " + genome + " sex: " + genomeRelationships.getSex(genome) + " output: " + outputGenomes.contains(genome));
        individualFactories.add(new IndividualSampleFactory<>(newParams, chooser, haploid, diploid, none, genomeRelationships.getSex(genome), sexMemo));
      }

      final PopulationCallerConfiguration pc = new PopulationCallerConfiguration(popCaller, outputGenomes.toArray(new String[0]), individualFactories, chooser, haploid, diploid, ssp, new ChildFamilyLookup(calledGenomes.size(), famArray));
      pc.getVcfAnnotators().addAll(annot);
      return pc;
    }

  }

  private final MendelianDenovoChecker mDenovoCorrect;

  private PopulationCallerConfiguration(MultisampleJointCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp, ChildFamilyLookup families) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
    mDenovoCorrect = new MendelianDenovoChecker(families);
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP, VcfFormatField.SCONT));
    return f;
  }

  @Override
  public DenovoChecker getDenovoChecker() {
    return mDenovoCorrect;
  }
}
