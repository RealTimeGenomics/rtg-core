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

package com.rtg.variant.bayes.multisample.statistics;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.reference.SexMemo;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.vcf.VariantStatistics;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;
import com.rtg.variant.bayes.multisample.IndividualSampleFactory;
import com.rtg.variant.bayes.multisample.JointCallerConfigurator;
import com.rtg.variant.bayes.multisample.MultisampleUtils;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;

/**
 */
public final class VstatsCallerConfiguration extends AbstractJointCallerConfiguration {

  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<VariantStatistics> {

    @Override
    public VstatsCallerConfiguration getConfig(final VariantParams params, VariantStatistics statistics) throws IOException {
      Diagnostic.userLog("Using singleton caller");

      final String[] outputSampleNames = SamUtils.getSampleNames(params.uberHeader());
      if (outputSampleNames.length > 1) {
        throw new NoTalkbackSlimException("Multiple samples detected in read group headers, but this command models only one genome");
      }

      final String sampleName = (outputSampleNames.length == 0) ? VariantOutputVcfFormatter.DEFAULT_SAMPLE : outputSampleNames[0];

      final String[] genomeNames = {sampleName};
      final VstatsCaller statisticsCaller = new VstatsCaller(params);
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
      individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, params.sex(), sexMemo));
      return new VstatsCallerConfiguration(statisticsCaller, genomeNames, individualFactories, chooser, haploid, diploid, ssp);
    }
  }

  VstatsCallerConfiguration(VstatsCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
  }

  /**
   * @throws IOException  when closing
   */
  @Override
  public void close() throws IOException {
    getJointCaller().close();
  }
}
