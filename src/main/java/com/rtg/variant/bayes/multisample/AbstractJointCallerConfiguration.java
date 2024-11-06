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

package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.DenovoChecker;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.multisample.population.SiteSpecificPriors;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfFilter;

/**
 * Warning: it is assumed that this class and its subclasses are immutable
 * and thus thread safe.
 */
@TestClass("com.rtg.variant.bayes.multisample.cancer.SomaticCallerConfigurationTest")
public abstract class AbstractJointCallerConfiguration {


  private final MultisampleJointCaller mJointCaller; // Makes the variant calls
  private final ModelSnpFactory mHaploid;            // For creation of hypotheses priors
  private final ModelSnpFactory mDiploid;            // For creation of hypotheses priors
  private final String[] mGenomeNames;
  private final List<IndividualSampleFactory<?>> mIndividualFactories;
  private final MachineErrorChooserInterface mMachineErrorChooser;
  private final PopulationHwHypothesesCreator mSiteSpecificPriors;
  private final Collection<VcfAnnotator> mAnnotators = new ArrayList<>();
  private final Collection<VcfFilter> mFilters = new ArrayList<>();

  protected AbstractJointCallerConfiguration(
      MultisampleJointCaller jointCaller,
      String[] genomeNames,
      List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser,
      ModelSnpFactory haploidFactory,
      ModelSnpFactory diploidFactory,
      PopulationHwHypothesesCreator siteSpecificPriors
  ) {
    mJointCaller = jointCaller;
    mGenomeNames = genomeNames;
    mIndividualFactories = individualFactories;
    mMachineErrorChooser = machineErrorChooser;
    mHaploid = haploidFactory;
    mDiploid = diploidFactory;
    mSiteSpecificPriors = siteSpecificPriors;
  }

  /**
   * @return site specific priors
   */
  public SiteSpecificPriors getSiteSpecificPriors() {
    return mSiteSpecificPriors;
  }

  /**
   * Gets names of all genomes for which calls will be output. There may be additional
   * genomes not included in this list (e.g. those present in a pedigree but for which
   * mappings were not provided).
   * @return names of all the genomes for which calls will be output.
   */
  public String[] getGenomeNames() {
    return mGenomeNames;
  }

  /**
   * @return machine error chooser
   */
  public MachineErrorChooserInterface getMachineErrorChooser() {
    return mMachineErrorChooser;
  }

  /**
   * Number of genomes for which calls are calculated. This may be higher than the number
   * of genomes that have calls output.
   * @return number of genomes.
   */
  public int numberOfGenomes() {
    return mIndividualFactories.size();
  }

  /**
   * Get output formatter
   * @param params input parameters
   * @return a new output formatter
   */
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    return new VariantOutputVcfFormatter(params, getGenomeNames());
  }

  /**
   * Test if this caller supports given ploidy
   * @param ploidy ploidy to test
   * @return true if this caller supports running on given ploidy
   */
  public boolean handlesPloidy(final Ploidy ploidy) {
    return true;
  }

  /**
   * Get the family (if any) in which given sample is a child.
   * @return the lookup table
   */
  public DenovoChecker getDenovoChecker() {
    return null;
  }

  /**
   * Get a list of any VcfAnnotators to apply to the output
   * @return a list of VcfAnnotators
   */
  public Collection<VcfAnnotator> getVcfAnnotators() {
    return mAnnotators;
  }

  /**
   * Get a list of any VcfFilters to apply to the output
   * @return a list of VcfFilters
   */
  public Collection<VcfFilter> getVcfFilters() {
    return mFilters;
  }

  // Ploidy can vary throughout the chromosome e.g. PAR regions
  protected Ploidy getEffectivePloidy(int sampleNumber, String refName, int pos) {
    return mIndividualFactories.get(sampleNumber).getEffectivePloidy(refName, pos);
  }

  /**
   * @return get joint caller
   */
  public MultisampleJointCaller getJointCaller() {
    return mJointCaller;
  }

  /**
   * Get a caller suitable for use with supplied hypotheses.
   *
   * @param complex hypotheses
   * @param params command line parameters and related information.
   * @param cot describes region on template that hypotheses replace.
   * @return joint caller
   */
  public MultisampleJointCaller getComplexJointCaller(Hypotheses<DescriptionComplex> complex, VariantParams params, ComplexTemplate cot) {
    return mJointCaller;
  }

  /**
   *
   * @param hyp hypotheses containing haploid and diploid priors
   * @param locus position to which the model should correspond
   * @return set of models for complex calling, one per sample.
   */
  public List<ModelInterface<?>> getModelComplex(HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> hyp, SequenceNameLocus locus) {
    final List<ModelInterface<?>> list = new ArrayList<>(mIndividualFactories.size());
    for (IndividualSampleFactory<?> individualFactory : mIndividualFactories) {
      list.add(individualFactory.makeModelComplex(hyp, locus));
    }
    return list;
  }

  /**
   * Gets sample processors for every sample.
   *
   * @param name of the sequence.
   * @param template nucleotides.
   * @param start position on the template of region being processed (0 based, inclusive).
   * @param end  position on the template of region being processed (0 based, exclusive).
   * @return a processor for each sample.
   */
  public IndividualSampleProcessor<?>[] getIndividualSampleProcessors(final String name, final byte[] template, int start, int end) {
    final IndividualSampleProcessor<?>[] res = new IndividualSampleProcessor<?>[mIndividualFactories.size()];
    for (int i = 0; i < res.length; ++i) {
      res[i] = mIndividualFactories.get(i).make(name, template, start, end);
    }
    return res;
  }

  /**
   * Gets the set of hypotheses for simple SNP calls with appropriate initial priors.
   * @param ref reference code id
   * @param templateName template name
   * @param pos position
   * @return the set of haploid and diploid hypotheses you'd use for the specified reference position
   */
  public HaploidDiploidHypotheses<HypothesesPrior<Description>> getSnpHypotheses(int ref, String templateName, int pos) {
    if (mSiteSpecificPriors != null) {
      final HaploidDiploidHypotheses<HypothesesPrior<Description>> ssp = mSiteSpecificPriors.getSnpHypotheses(templateName, pos);
      if (ssp != null) {
        return ssp;
      }
//      if (PopulationHwHypothesesCreator.NOVELTY_FACTOR > 0) {
//        return PopulationHwHypothesesCreator.computeNoveltyPriors(new HaploidDiploidHypotheses<Hypotheses<?>>(HypothesesNone.SINGLETON, mHaploid.make(ref, templateName, pos).hypotheses(), mDiploid.make(ref, templateName, pos).hypotheses()), PopulationHwHypothesesCreator.NOVELTY_FACTOR);
//      }
    }
    return new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, mHaploid.defaultHypotheses(ref), mDiploid.defaultHypotheses(ref));
  }

  /**
   * Called at end of all processing.
   * @throws IOException whenever.
   */
  public void close()  throws IOException {
    //by default do nothing.
  }
}
