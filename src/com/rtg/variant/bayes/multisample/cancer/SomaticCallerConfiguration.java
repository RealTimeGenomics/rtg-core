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
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.reference.Ploidy;
import com.rtg.reference.SexMemo;
import com.rtg.relation.LineageLookup;
import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.sam.SamUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.SomaticParams;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.DenovoChecker;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.complex.LineageDenovoChecker;
import com.rtg.variant.bayes.complex.StatisticsComplex;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.IndividualSampleFactory;
import com.rtg.variant.bayes.multisample.JointCallerConfigurator;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.MultisampleUtils;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;
import com.rtg.variant.realign.AlignmentEnvironment;
import com.rtg.variant.realign.AlignmentEnvironmentGenomeSubstitution;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.EnvironmentCombined;
import com.rtg.variant.realign.RealignParamsGenome;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 */
public final class SomaticCallerConfiguration extends AbstractJointCallerConfiguration {


  /**
   * The factory for this caller.
   */
  public static final class Configurator implements JointCallerConfigurator<SomaticStatistics> {

    @Override
    public SomaticCallerConfiguration getConfig(final VariantParams params, final SomaticStatistics statistics) throws IOException {
      final double loh = params.somaticParams().lohPrior();
      final int numberOfGenomes = params.genomeRelationships().genomes().length;
      final Relationship[] derived = params.genomeRelationships().relationships(RelationshipType.ORIGINAL_DERIVED);
      assert derived.length == 1 || numberOfGenomes == 2;
      final String[] genomeNames = {derived[0].first(), derived[0].second()};
      final String normalSampleName = genomeNames[AbstractSomaticCaller.NORMAL];
      final String cancerSampleName = genomeNames[AbstractSomaticCaller.CANCER];
      final String[] outputSampleNames = SamUtils.getSampleNames(params.uberHeader());
      if (outputSampleNames.length != 2) {
        throw new NoTalkbackSlimException("Exactly two sample names expected in mappings");
      }
      for (final String mapName : outputSampleNames) {
        if (!mapName.equals(normalSampleName) && !mapName.equals(cancerSampleName)) {
          throw new NoTalkbackSlimException("Unexpected sample name in mappings: " + mapName);
        }
      }
      final Double contam = derived[0].getContamination();
      if (contam == null) {
        throw new RuntimeException("Contamination level not specified in genome relationship file or on command line.");
      }
      final double contamination = contam;
      final Double reverseContam = derived[0].getReverseContamination();
      final double reverseContamination = reverseContam == null ? 0 : reverseContam;
      final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
      final PopulationHwHypothesesCreator ssp;
      final double contraryProbability = params.genomePriors().contraryProbability();
      final double phi = Math.min(1, contraryProbability + reverseContamination);
      final double psi = Math.min(1, contraryProbability + contamination);
      if (params.populationPriorFile() != null) {
        ssp = new PopulationHwHypothesesCreator(params.populationPriorFile(), params.genomePriors(), params.referenceRanges(), params.alleleBalance());
      } else {
        ssp = null;
      }
      final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, params.alleleBalance());
      final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, params.alleleBalance());
      final ModelNoneFactory none = new ModelNoneFactory();
      final List<IndividualSampleFactory<?>> individualFactories = new ArrayList<>();
      final SexMemo sexMemo = Utils.createSexMemo(params);
      final AbstractSomaticCaller jointCaller;
      if (contamination == 0.0) {
        Diagnostic.userLog("Using no contamination cancer caller");
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, params.sex(), sexMemo));
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, params.sex(), sexMemo));
        jointCaller = new PureSomaticCaller(
          new CachedSomaticPriorsFactory<>(haploid.defaultHypotheses(0), loh),
          new CachedSomaticPriorsFactory<>(diploid.defaultHypotheses(0), loh),
          params, phi, psi);
      } else {
        Diagnostic.userLog("Using contaminated cancer caller");
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, params.sex(), sexMemo));
        final ModelCancerFactory contamHaploid = new ModelCancerFactory(params.genomePriors(), contamination, true, params.alleleBalance());
        final ModelCancerFactory contamDiploid = new ModelCancerFactory(params.genomePriors(), contamination, false, params.alleleBalance());
        individualFactories.add(new IndividualSampleFactory<>(params, chooser, contamHaploid, contamDiploid, none, params.sex(), sexMemo));
        jointCaller = new ContaminatedSomaticCaller(
          new CachedSomaticPriorsFactory<>(haploid.defaultHypotheses(0), loh),
          new CachedSomaticPriorsFactory<>(diploid.defaultHypotheses(0), loh),
          params, phi, psi);
      }
      final SomaticCallerConfiguration sc = new SomaticCallerConfiguration(jointCaller, genomeNames, individualFactories, chooser, contamination, haploid, diploid, ssp, phi, psi, params.alleleBalance());
      sc.getVcfFilters().add(new SomaticFilter(statistics, !(params.somaticParams().includeGermlineVariants() || params.callLevel() == VariantOutputLevel.ALL)));
      return sc;
    }
  }

  private final double mContamination;
  private final double mPhi;
  private final double mPsi;
  private final LineageDenovoChecker mChecker;
  private final AlleleBalanceProbability mAlleleBalanceProbability;

  private SomaticCallerConfiguration(MultisampleJointCaller jointCaller, String[] genomeNames, List<IndividualSampleFactory<?>> individualFactories, MachineErrorChooserInterface machineErrorChooser, double contamination, ModelSnpFactory haploid, ModelSnpFactory diploid, PopulationHwHypothesesCreator ssp, double phi, double psi, AlleleBalanceProbability alleleProbability) {
    super(jointCaller, genomeNames, individualFactories, machineErrorChooser, haploid, diploid, ssp);
    mContamination = contamination;
    mPhi = phi;
    mPsi = psi;
    mChecker = new LineageDenovoChecker(new LineageLookup(-1, 0)); // normal is sample 0, cancer is sample 1
    mAlleleBalanceProbability  = alleleProbability;
  }

  @Override
  public DenovoChecker getDenovoChecker() {
    return mChecker;
  }

  @Override
  public VariantOutputVcfFormatter getOutputFormatter(final VariantParams params) {
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(params, getGenomeNames());
    f.addExtraInfoFields(EnumSet.of(VcfInfoField.LOH, VcfInfoField.NCS));
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.SSC, VcfFormatField.SS, VcfFormatField.COC, VcfFormatField.COF, VcfFormatField.VAF));
    return f;
  }

  @Override
  public List<ModelInterface<?>> getModelComplex(HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> hyp, SequenceNameLocus locus) {
    if (mContamination == 0.0) {
      return super.getModelComplex(hyp, locus);
    }
    assert numberOfGenomes() == 2;
    final List<ModelInterface<?>> list = new ArrayList<>();
    final Ploidy ploid = getEffectivePloidy(0, locus.getSequenceName(), (locus.getStart() + locus.getEnd()) / 2);
    final Hypotheses<DescriptionComplex> complex;
    switch (ploid) {
      case HAPLOID:
        complex = hyp.haploid();
        break;
      default:
        complex = hyp.diploid();
        break;
    }
    list.add(new Model<>(complex, new StatisticsComplex(complex.description(), locus.getLength()), mAlleleBalanceProbability));
    list.add(new ModelCancerContamination<>(new HypothesesCancer<>(complex, LogApproximatePossibility.SINGLETON), mContamination, new StatisticsComplex(complex.description(), locus.getLength()), mAlleleBalanceProbability));
    return list;
  }

  @Override
  public MultisampleJointCaller getComplexJointCaller(final Hypotheses<DescriptionComplex> complex, final VariantParams params, final ComplexTemplate cot) {
    final double[][] initialPriors = makeSomaticInitialPriors(complex.description(), cot);
    //System.err.println("CancerJointCallerConfiguration mu=" + mMutationRate + " contamination=" + mContamination);
    //System.err.println(complex);
    //System.err.println(toStringLog(initialPriors));
    // Currently this is just computing one Q corresponding to the hypotheses passed in -- it would be better to do both if we want to do multiple ploidys at once in a multiple model situation.
    final SomaticParams somaticParams = params.somaticParams();
    final SomaticPriorsFactory<DescriptionComplex> qHaploidFactory = complex.haploid() ? new SomaticPriorsFactory<>(complex, somaticParams.lohPrior(), initialPriors) : null;
    final SomaticPriorsFactory<DescriptionComplex> qDiploidFactory = complex.haploid() ? null : new SomaticPriorsFactory<>(complex, somaticParams.lohPrior(), initialPriors);
    //System.err.println(toStringLog(q));
    if (mContamination == 0) {
      return new PureSomaticCaller(qHaploidFactory, qDiploidFactory, params, mPhi, mPsi);
    } else {
      return new ContaminatedSomaticCaller(qHaploidFactory, qDiploidFactory, params, mPhi, mPsi);
    }
  }

  /**
   * This diagram shows the relationship of the different variables when doing the all-paths matching.
   * <img src="doc-files/makeInitialPriors.jpg" alt="image">
   * @param description haploid hypotheses.
   * @param cot region on reference being replaced by the hypotheses.
   * @return the prior probabilities for transitions between haploid hypotheses.
   */
  static double[][] makeSomaticInitialPriors(final Description description, final ComplexTemplate cot) {
    final int hypExtension = Math.max(5, Math.max(description.maxLength() + 1, cot.getLength() + 1));
    final AllPaths sm = new ScoreFastUnderflow(RealignParamsGenome.SINGLETON);
    final double[][] initialPriorsLn = new double[description.size()][description.size()];
    final int start = cot.getStart() - hypExtension;
    final int end = cot.getStart() + hypExtension;
    //System.err.println("CancerJointCallerConfiguration start=" + start + " end=" + end);
    final int e0Hx = cot.getEnd() - hypExtension;
    final byte[][] callBytes = new byte[description.size()][];
    for (int i = 0; i < description.size(); i++) {
      callBytes[i] = DNA.stringDNAtoByte(description.name(i));
    }
    for (int i = 0; i < description.size(); i++) {
      final int iLen = callBytes[i].length;
      final AlignmentEnvironment aei = new AlignmentEnvironmentGenomeSubstitution(start, 0 /*doesn't matter*/, cot, callBytes[i]);
      //System.err.println("aei.length=" + aei.length());
      for (int j = 0; j < description.size(); j++) {
        final int jLen = callBytes[j].length;
        final AlignmentEnvironment aej = new AlignmentEnvironmentGenomeSubstitution(start, end, cot, callBytes[j]);
        final int s = e0Hx - (iLen + jLen) / 2;
        final int mx = (iLen + jLen + 1) / 2;
        final EnvironmentCombined env = new EnvironmentCombined(aej, s, mx, aei);
        sm.setEnv(env);
        initialPriorsLn[i][j] = sm.totalScoreLn();
        //System.err.println(sm.toString());
      }
    }
    //System.err.println(IntegralAbstract.toString(initialPriorsLn));
    final double[][] initialPriors = new double[description.size()][];
    for (int i = 0; i < initialPriorsLn.length; i++) {
      initialPriors[i] = MathUtils.lnToNormaliedProb(initialPriorsLn[i]);
    }
    return initialPriors;
  }

}
