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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.rtg.launcher.GlobalFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.reference.Ploidy;
import com.rtg.sam.ReaderWindow;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.complex.EvidenceComplex;
import com.rtg.variant.bayes.complex.HypothesesComplex;
import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;
import com.rtg.variant.bayes.multisample.population.DescriptionCounts;
import com.rtg.variant.bayes.multisample.population.HwEstimator;
import com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreator;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 */
final class ComplexCaller {

  // The maximum number of hypotheses that can comfortably be handled by the complex caller.
  // Calling will bail out where this is exceeded.
  // Usually the set of hypotheses is sufficiently pruned prior to this point that
  // the cutoff only comes into play to prevent very rare pathological cases.
  private static final int MAX_HYPOTH = GlobalFlags.getIntegerValue(GlobalFlags.COMPLEX_CALLER_MAX_HYPOTH_FLAG);

  private static final int MAX_HYPOTHESES = (int) Math.sqrt(Integer.MAX_VALUE);

  private int mExcessCoverageCount = 0;
  private int mNoHypothesesCount = 0;
  private int mExcessHypothesesCount = 0;

  private final VariantParams mParams;
  private final AbstractJointCallerConfiguration mConfig;

  ComplexCaller(VariantParams params,  AbstractJointCallerConfiguration config) {
    mParams = params;
    mConfig = config;
  }

  List<Variant> makeComplexCalls(Complexities complexities, ReaderWindow<VariantAlignmentRecord> tribble, byte[] template, String refName)  throws IOException {
    assert complexities.isFixed();
    final ArrayList<Variant> calls = new ArrayList<>();
    final int numSamples = mConfig.numberOfGenomes();

    for (final ComplexRegion region : complexities) {
      if (region.type() == RegionType.HYPER || region.type() == RegionType.OVERFLOW) {
        //bed file will record this
        continue;
      }
      final int startOfRegion = region.getStart();
      final int endOfRegion = region.getEnd();

      // Need to check the region is within the bounds of the specified complexities
      // and if the start is within the complex region then only process it
      if (startOfRegion >= complexities.startOfChunk() && startOfRegion <= complexities.endOfChunk()) {
        final ComplexTemplate cot = new ComplexTemplate(template, refName, startOfRegion, endOfRegion);

        final Iterator<VariantAlignmentRecord> it = tribble.recordsOverlap(startOfRegion, endOfRegion);
        final List<AlignmentMatch> matches = MultisampleUtils.intersectSet(startOfRegion, endOfRegion, it, mConfig.getMachineErrorChooser(), mParams);

        // Create the complex hypotheses
        final HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> hyp = makeComplexHypotheses(region, cot, matches, calls);
        if (hyp == null) { // Couldn't create suitable hypotheses, makeComplexHypotheses will have updated calls and region appropriately
          continue;
        }

        // Pass evidence to correct models
        final List<ModelInterface<?>> models = mConfig.getModelComplex(hyp, cot);
        for (final AlignmentMatch match : matches) {
          final int genome = (numSamples == 1) ? 0 : match.alignmentRecord().getGenome();
          if (!(models.get(genome) instanceof ModelNone)) {
            final HypothesesPrior<DescriptionComplex> hypotheses = models.get(genome).haploid() ? hyp.haploid() : hyp.diploid();
            final EvidenceComplex evidence = new EvidenceComplex(hypotheses, match, cot, mParams, mConfig.getMachineErrorChooser());
            models.get(genome).increment(evidence);
          }
        }


        // Make the calls
        final MultisampleJointCaller complexJointCaller = mConfig.getComplexJointCaller(hyp.get(models.get(0)), mParams, cot); /// The first argument is only needed by the cancer caller -- seems bogus.
        final Variant variant;
        for (ModelInterface<?> model : models) {
          model.freeze();
        }
        final Variant variant0 = complexJointCaller.makeCall(refName, startOfRegion, endOfRegion, template, models, hyp);

        if (variant0 == null && mParams.callLevel() == VariantOutputLevel.ALL) {
          //need to output an identity call at this position
          final VariantLocus locus = new VariantLocus(refName, startOfRegion, endOfRegion, DnaUtils.bytesToSequenceIncCG(template, startOfRegion, endOfRegion - startOfRegion), VariantUtils.getPreviousRefNt(template, startOfRegion));
          final VariantSample[] samples = new VariantSample[numSamples];
          for (int i = 0; i < numSamples; i++) {
            samples[i] = new VariantSample(Ploidy.NONE);
          }
          variant = new Variant(locus, samples);
        } else {
          variant = variant0;
        }
        if (variant != null) {
          variant.setComplexScored();
          if (variant.hasCallNames()) {
            calls.add(variant);
          } else {
            copyFailedCalls(complexities, calls, startOfRegion, endOfRegion, mParams); //mark the calls as hyper complex and copy over the originals in complex call list
          }
        } else {
          region.setType(RegionType.COMPLEX_NO_VARIANT);
        }
      }
    }
    return calls;
  }

  private HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> makeComplexHypotheses(ComplexRegion region, ComplexTemplate cot, List<AlignmentMatch> matches, ArrayList<Variant> calls) {
    final String refName = cot.getSequenceName();
    final int overcoverageShortCircuit = mParams.maxCoverageBypass().thresholdTotal(refName);
    final HwEstimator hwEstimator = new HwEstimator();

    if (matches == null) { // overflow encountered
      // complex SSP is forcing complex calling in over coverage region, mark it as over coverage and continue to next one
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, VariantFilter.COVERAGE.mask(), mParams.maxCoverageBypass().thresholdTotal(refName)));
      }
      mExcessCoverageCount++;
      region.setType(RegionType.OVERCOVERAGE);
      return null;
    }

    //System.err.println("CC: [" + startOfRegion + "-" + endOfRegion + ") " + matches.size() + " overlapping matches");
    //for (AlignmentMatch match : matches) {
    //  System.err.println("Match [" + match.alignmentRecord().getStart() + "-" + match.alignmentRecord().getEnd() + ") " + match.alignmentRecord().getCigar());
    //}

    // If the coverage is too high make a dummy call that indicates over-coverage situation
    if (matches.size() > overcoverageShortCircuit) {
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, VariantFilter.COVERAGE.mask(), matches.size()));
      }
      mExcessCoverageCount++;
      region.setType(RegionType.OVERCOVERAGE);
      return null;
    }
    if (!canCreateHypothesis(matches)) {
      Diagnostic.userLog("Can not create hypothesis to cover region " + cot.toString());
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, 0, matches.size()));
      }
      mNoHypothesesCount++;
      region.setType(RegionType.NO_HYPOTHESES);
      return null;
    }

    cot.setComplexContext(HypothesesComplex.createComplexDescription(matches, cot, mConfig.getSiteSpecificPriors(), mParams.pruneHypotheses(), mParams.maxComplexHypotheses()), LogApproximatePossibility.SINGLETON);

    final HypothesesComplex hypHap = HypothesesComplex.makeComplexHypotheses(cot, true, mParams);
    if (hypHap.size() == 0 || (hypHap.size() == 1 && hypHap.reference() == 0)) {
      Diagnostic.userLog("Can not create hypothesis to cover region " + cot.toString());
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, 0, matches.size()));
      }
      mNoHypothesesCount++;
      region.setType(RegionType.NO_HYPOTHESES);
      return null;
    }
    if (mParams.pruneHypotheses() && hypHap.size() > MAX_HYPOTH) { // This is after pruning, so will almost never actually trigger.
      Diagnostic.userLog(hypHap.size() + " haploid hypotheses is too high, at region " + cot.toString());
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, 0, matches.size()));
      }
      mExcessHypothesesCount++;
      region.setType(RegionType.TOO_MANY_HYPOTHESES);
      return null;
    }

    final HypothesesComplex hypDip = HypothesesComplex.makeComplexHypotheses(cot, false, mParams);
    if (hypHap.size() > MAX_HYPOTHESES || hypDip.size() > MAX_HYPOTHESES) {
      Diagnostic.userLog(hypHap.size() + " haploid hypotheses, " + hypDip.size() + " diploid hypotheses is too high, at region " + cot.toString());
      if (mParams.callLevel() == VariantOutputLevel.ALL) {
        calls.add(makeDummyCall(cot, 0, matches.size()));
      }
      mExcessHypothesesCount++;
      region.setType(RegionType.TOO_MANY_HYPOTHESES);
      return null;
    }
    final DescriptionCounts dc = mConfig.getSiteSpecificPriors() == null ? null : PopulationHwHypothesesCreator.getDescriptionCounts(hypHap.description(), hypHap.reference());
    return PopulationHwHypothesesCreator.createHapDipHypotheses(HypothesesNone.SINGLETON_COMPLEX, hypHap, hypDip, dc, hwEstimator);
  }

  private void copyFailedCalls(Complexities complexities, final ArrayList<Variant> calls, final int startOfRegion, final int endOfRegion, VariantParams params) {
    for (final Variant v : complexities.getOriginalCalls()) {
      if (v != null && v.getLocus().getStart() >= startOfRegion && v.getLocus().getEnd() <= endOfRegion) {
        if (params.callLevel() == VariantOutputLevel.ALL && !v.isIndel() || (v.getNumberOfSamples() > 0 && v.getSample(0) != null && v.getSample(0).getName() != null && !v.getSample(0).isIdentity())) {
          v.addFilter(VariantFilter.HYPER_COMPLEX);
          calls.add(v);
        }
      }
    }
  }

  private static boolean canCreateHypothesis(List<AlignmentMatch> matches) {
    // TODO Take into account any SSP hypotheses?
    boolean res = false;
    for (final AlignmentMatch m : matches) {
      if (m.isFixedLeft() && m.isFixedRight()) {
        res = true;
      }
    }
    return res;
  }

  // make a dummy call that indicates a failure situation
  private Variant makeDummyCall(ComplexTemplate cot, int filterMask, int numMatches) {
    final int numSamples = mConfig.numberOfGenomes();
    final VariantSample[] samples = new VariantSample[numSamples];
    for (int k = 0; k < numSamples; k++) {
      //SAI: matches.size() here is possibly wrong for multisample cases, but
      //I don't think it is very important, we only output these with --all mode
      //and these aren't really proper "calls" anyway.
      final VariantSample vs = new VariantSample(mConfig.getEffectivePloidy(k, cot.getSequenceName(), (cot.getStart() + cot.getEnd()) / 2));
      vs.setCoverage(numMatches, 0.0);
      samples[k] = vs;
    }
    final VariantLocus locus = new VariantLocus(cot.getSequenceName(), cot.getStart(), cot.getEnd(), cot.replaceString(), VariantUtils.getPreviousRefNt(cot.templateBytes(), cot.getStart()));
    final Variant v = new Variant(locus, filterMask, samples);
    v.setComplexScored();
    v.addFilter(VariantFilter.FAILED_COMPLEX);
    return v;
  }

  /**
   * @return number of calls that were not made due to excessively high coverage.
   */
  public int getExcessiveCoverageCount() {
    return mExcessCoverageCount;
  }

  /**
   * @return number of calls that were not made due to excessively number of hypotheses.
   */
  public int getExcessiveHypothesesCount() {
    return mExcessHypothesesCount;
  }

  /**
   * @return number of calls that were not made due to having no hypotheses (probably due to no reads spanning region).
   */
  public int getNoHypothesesCount() {
    return mNoHypothesesCount;
  }
}
