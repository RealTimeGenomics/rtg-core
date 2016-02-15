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
package com.rtg.variant.bayes.multisample;


import com.rtg.reference.Ploidy;
import com.rtg.reference.Sex;
import com.rtg.reference.SexMemo;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.SamToMatch;
import com.rtg.variant.SamToMatchCigar;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.complex.StatisticsComplex;
import com.rtg.variant.bayes.snp.CigarParserModel;
import com.rtg.variant.bayes.snp.EvidenceMatcher;
import com.rtg.variant.bayes.snp.EvidenceQFactory;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.IndelDetector;
import com.rtg.variant.bayes.snp.IndelDetectorFactory;
import com.rtg.variant.bayes.snp.IndelIgnorerFactory;
import com.rtg.variant.bayes.snp.IndelMatcher;
import com.rtg.variant.bayes.snp.ReferenceBasedBuffer;
import com.rtg.variant.bayes.snp.SwitchingReferenceBasedBuffer;

/**
 * A wrapper class around a <code>BayesMatcher</code> object.
 * This class chooses the appropriate kind of individual caller
 * for the current reference template.
 * @param <D> description type that these will generate for haploid and diploid callers
 */
public class IndividualSampleFactory<D extends Description> {

  // If true, call appropriate ploidy in PAR regions.
  static final boolean PAR_AWARE = true; //GlobalFlags.getBooleanValue(GlobalFlags.CALLER_PAR_AWARE);

  private final VariantParams mParams;
  private final MachineErrorChooserInterface mChooser;
  private final ModelFactory<D, ?> mHaploidFactory;
  private final ModelFactory<D, ?> mDiploidFactory;
  private final ModelFactory<D, ?> mNoneFactory;
  private final Sex mSex;
  private final SexMemo mSexMemo;

  /**
   * Construct a new variant task with given params and output.
   * @param params parameters
   * @param chooser used to select machine errors for a SAM record.
   * @param haploidFactory factory used for haploid sequences.
   * @param diploidFactory factory used for diploid sequences.
   * @param noneFactory factory used for none sequences.
   * @param sex of the sample.
   * @param sexMemo information from reference file and parameters that determines ploidy given sequence and sex.
   */
  public IndividualSampleFactory(final VariantParams params, final MachineErrorChooserInterface chooser, final ModelFactory<D, ?> haploidFactory, final ModelFactory<D, ?> diploidFactory, final ModelFactory<D, ?> noneFactory, final Sex sex, SexMemo sexMemo) {
    mParams = params;
    mChooser = chooser;
    mHaploidFactory = haploidFactory;
    mDiploidFactory = diploidFactory;
    mNoneFactory = noneFactory;
    mSex = sex;
    mSexMemo = sexMemo;
  }

  Ploidy getEffectivePloidy(String refName, int pos) {
    return PAR_AWARE ? mSexMemo.getEffectivePloidy(mSex, refName, pos) : mSexMemo.getEffectivePloidy(mSex, refName);
  }

  // Used by complex calling
  ModelInterface<DescriptionComplex> makeModelComplex(HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> hyp, SequenceNameLocus locus) {
    switch (getEffectivePloidy(locus.getSequenceName(), (locus.getStart() + locus.getEnd()) / 2)) {
      case HAPLOID:
        return new Model<>(hyp.haploid(), new StatisticsComplex(hyp.haploid().description(), locus.getLength()), mParams.alleleBalance());
      case NONE:
        return ModelNone.SINGLETON_COMPLEX;
      default:
        return new Model<>(hyp.diploid(), new StatisticsComplex(hyp.haploid().description(), locus.getLength()), mParams.alleleBalance());
    }
  }

  // Used by simple snp calling
  IndividualSampleProcessor<?> make(final String name, final byte[] template, int start, int end) {
    final ReferenceBasedBuffer<ModelInterface<D>> result;
    final ReferenceBasedBuffer<IndelDetector> indelDetector;
    final SequenceNameLocusSimple region = new SequenceNameLocusSimple(name, start, end);
    final int parBoundary = PAR_AWARE ? mSexMemo.getParBoundary(mSex, region) : -1;
    if (parBoundary == -1) {
      result = new ReferenceBasedBuffer<>(end - start, selectFactory(name, start), template, start);
      indelDetector = new ReferenceBasedBuffer<>(end - start, IndelDetectorFactory.SINGLETON, template, start);
    } else {
      // Make a special buffer that switches factories when it crosses the PAR boundary. The assumption is that there is only ever one such boundary within a chunk.
      Diagnostic.developerLog("Creating a " + mSex + " PAR-aware boundary chunk " + region + " crossing at " + (parBoundary + 1)
          + " from " + mSexMemo.getEffectivePloidy(mSex, name, start)
          + " to " + mSexMemo.getEffectivePloidy(mSex, name, end));
      result = new SwitchingReferenceBasedBuffer<>(end - start, selectFactory(name, start), selectFactory(name, end), parBoundary, template, start);
      indelDetector = new SwitchingReferenceBasedBuffer<>(end - start, selectIndelFactory(name, start), selectIndelFactory(name, end), parBoundary, template, start);
    }
    final EvidenceMatcher<ModelInterface<D>> matcherCurrent = new EvidenceMatcher<>(result, new EvidenceQFactory());
    final IndelMatcher indelMatcher = new IndelMatcher(indelDetector);
    final SamToMatch toMatch = new SamToMatchCigar(mParams, new CigarParserModel(matcherCurrent, indelMatcher, start, end, mParams), mChooser);
    return new IndividualSampleProcessor<>(template, matcherCurrent, indelMatcher, toMatch);
  }

  private ModelFactory<D, ?> selectFactory(final String name, int pos) {
    switch (getEffectivePloidy(name, pos)) {
      case HAPLOID:
        return mHaploidFactory;
      case NONE:
        return mNoneFactory;
      default:
        return mDiploidFactory;
    }
  }

  private IndelDetectorFactory selectIndelFactory(final String name, int pos) {
    switch (getEffectivePloidy(name, pos)) {
      case HAPLOID:
        return IndelDetectorFactory.SINGLETON;
      case NONE:
        return IndelIgnorerFactory.SINGLETON;
      default:
        return IndelDetectorFactory.SINGLETON;
    }
  }
}
