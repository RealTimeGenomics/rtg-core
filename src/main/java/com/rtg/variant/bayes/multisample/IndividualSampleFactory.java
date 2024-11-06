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
