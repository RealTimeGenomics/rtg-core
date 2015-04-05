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
package com.rtg.variant.eval;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IORunnable;
import com.rtg.util.Pair;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;

/**
 * Runs all the evaluation for a single reference sequence
 */
@TestClass({"com.rtg.variant.eval.VcfEvalTaskTest", "com.rtg.variant.eval.SequenceEvaluatorTest"})
class SequenceEvaluator implements IORunnable {

  private final EvalSynchronizer mSynchronize;
  private final SequencesReader mTemplate;
  private final Map<String, Long> mNameMap;

  SequenceEvaluator(EvalSynchronizer variantSets, Map<String, Long> nameMap, SequencesReader template) {
    mSynchronize = variantSets;
    mTemplate = template;
    mNameMap = nameMap;
  }
  @Override
  public void run() throws IOException {
    final Pair<String, Map<VariantSetType, List<DetectedVariant>>> setPair = mSynchronize.nextSet();
    if (setPair == null) {
      return;
    }
    final String currentName = setPair.getA();
    final Long sequenceId = mNameMap.get(currentName);
    if (sequenceId == null) {
      throw new NoTalkbackSlimException("Sequence " + currentName + " is not contained in the reference.");
    }
    final byte[] template = mTemplate.read(sequenceId);

    final Map<VariantSetType, List<DetectedVariant>> set = setPair.getB();

    final Collection<DetectedVariant> baseLineCalls = set.get(VariantSetType.BASELINE);
    final Collection<DetectedVariant> calledCalls = set.get(VariantSetType.CALLS);

    if (baseLineCalls == null || baseLineCalls.size() == 0) {
      mSynchronize.write(currentName, null, calledCalls, null, null);
      if (calledCalls != null) {
        Diagnostic.developerLog("Number of called variants: " + calledCalls.size());

          for (final DetectedVariant v : calledCalls) {
            mSynchronize.addRocLine(new RocLine(v.getSequenceName(), v.getStart(), v.getSortValue(), 0.0, false), v);
          }

      }
    } else if (calledCalls == null || calledCalls.size() == 0) {
      Diagnostic.developerLog("Number of baseline variants: " + baseLineCalls.size());
      mSynchronize.write(currentName, null, null, baseLineCalls, null);
    } else {

      Diagnostic.developerLog("Sequence: " + currentName + " has " + baseLineCalls.size() + " baseline variants");
      Diagnostic.developerLog("Sequence: " + currentName + " has " + calledCalls.size() + " called variants");

      //find the best path for variant calls
      final Path best = PathFinder.bestPath(template, currentName, calledCalls, baseLineCalls);
      //System.out.println(path);
      List<OrientedVariant> truePositives = best.getCalledIncluded();
      final List<OrientedVariant> baselineTruePositives = best.getBaselineIncluded();

      final List<Variant> falsePositives = best.getCalledExcluded();
      final List<Variant> falseNegatives = best.getBaselineExcluded();

      mSynchronize.addVariants(baselineTruePositives.size() + falseNegatives.size());
      Diagnostic.developerLog("Writing variants...");
      final Pair<List<OrientedVariant>, List<OrientedVariant>> calls = Path.calculateWeights(best, truePositives, baselineTruePositives);
      final PhasingResult misPhasings = countMisphasings(best);
      mSynchronize.addPhasings(misPhasings.mMisPhasings, misPhasings.mCorrectPhasings, misPhasings.mUnphaseable);

      // this step is currently necessary as sometimes you can (rarely) have variants included in the best path
      // but they do not correspond to any variant in baseline.
      // E.g. two variants which when both replayed cancel each other out.
      truePositives = calls.getA();
      merge(falsePositives, calls.getB());

      mSynchronize.write(currentName, truePositives, falsePositives, falseNegatives, baselineTruePositives);
      Diagnostic.developerLog("Generating ROC data...");

      double tpTotal = 0.0;
      for (final OrientedVariant v : truePositives) {
        final DetectedVariant dv = (DetectedVariant) v.variant();
        mSynchronize.addRocLine(new RocLine(dv.getSequenceName(), dv.getStart(), dv.getSortValue(), v.getWeight(), true), dv);
        tpTotal += v.getWeight();
      }
      if (tpTotal - baselineTruePositives.size() > 0.001) {
        throw new SlimException("true positives does not match baseline number tp weighted= " + tpTotal + " baseline = " + baselineTruePositives.size());
      }
      for (final Variant v : falsePositives) {
        final DetectedVariant dv = (DetectedVariant) (v instanceof OrientedVariant ? ((OrientedVariant) v).variant() : v);
        // System.out.println(dv);
        mSynchronize.addRocLine(new RocLine(dv.getSequenceName(), dv.getStart(), dv.getSortValue(), 0.0, false), dv);
      }
    }
  }
  static boolean groupInPhase(List<VariantSummary> group) {
    if (group.size() < 1) {
      return true;
    }
    if (!group.get(0).isPhased()) {
      return false;
    }
    final boolean phase = group.get(0).phase();
    for (VariantSummary v : group) {
      if (!v.isPhased() || v.phase() != phase) {
        return false;
      }
    }
    return true;
  }

  static PhasingResult countMisphasings(Path best) {
    final List<OrientedVariant> baselineIncluded = best.getBaselineIncluded();
    final List<OrientedVariant> callsIncluded = best.getCalledIncluded();
    final List<Variant> baselineExcluded = best.getBaselineExcluded();
    final List<Variant> callsExcluded = best.getCalledExcluded();

    final CallIterator baseline = new CallIterator(baselineIncluded, baselineExcluded);
    final CallIterator calls = new CallIterator(callsIncluded, callsExcluded);
    final List<Path.SyncPoint> sync = best.getSyncPoints();

    int misPhasings = 0;
    int unphaseable = 0;
    int correctPhasings = 0;
    boolean baseIsPhased = false;
    boolean basePhase = false;
    boolean callIsPhased = false;
    boolean callPhase = false;
    VariantSummary call = null;
    VariantSummary base = null;
    // Rather than assuming baseline is all phased, we'll be a bit more careful.
    for (Path.SyncPoint point : sync) {
      final List<VariantSummary> baselineSection = new ArrayList<>();
      do {
        if (base != null && base.startPos() < point.getPos()) {
          baselineSection.add(base);
          base = baseline.hasNext() ? baseline.next() : null;
        }
        if (base == null) {
          base = baseline.hasNext() ? baseline.next() : null;
        }
      } while (base != null && base.startPos() < point.getPos());
      final List<VariantSummary> callSection = new ArrayList<>();
      do {
        if (call != null) {
          callSection.add(call);
        }
        call = calls.hasNext() ? calls.next() : null;
      } while (call != null && call.startPos() < point.getPos());

      // We need all of these to be in the same phasing otherwise we can't tell which calls are swapped
      if (!groupInPhase(baselineSection)) {
        for (VariantSummary summary : callSection) {
          if (summary.isPhased()) {
            unphaseable++;
          }
        }
        baseIsPhased = false;
        callIsPhased = false;
        continue;
      }
      // When the baseline calls have flipped orientation
      boolean transition = false;
      if (!baseIsPhased) {
        callIsPhased = false;
      }

      for (VariantSummary baseSummary : baselineSection) {
        if (baseSummary.isPhased()) {
          if (basePhase != baseSummary.phase()) {
            transition = true;
          }
          baseIsPhased = true;
        }
        basePhase = baseSummary.phase();
      }

      for (VariantSummary currentCall : callSection) {
        if (!currentCall.isPhased()) {
          callIsPhased = false;
        } else if (!callIsPhased) {
          //Start phasing run
          callIsPhased = true;
          callPhase = currentCall.phase();
        } else {
          //Continue phasing
          final boolean callTransition = currentCall.phase() != callPhase;
          if (currentCall.included() && !(callTransition == transition)) {
            misPhasings++;
            callPhase = currentCall.phase();
          } else if (currentCall.included()) {
            correctPhasings++;
            callPhase = currentCall.phase();
          }
        }
        transition = false;
      }
    }
    return new PhasingResult(misPhasings, correctPhasings, unphaseable);
  }
  static class PhasingResult {
    final int mMisPhasings;
    final int mCorrectPhasings;
    final int mUnphaseable;

    PhasingResult(int misPhasings, int correctPhasings, int unphaseable) {
      mMisPhasings = misPhasings;
      mCorrectPhasings = correctPhasings;
      mUnphaseable = unphaseable;
    }
  }

  /**
   * Combines included/excluded calls into one stream for purposes of bridging phasing across fp
   */
  private static class CallIterator implements Iterator<VariantSummary> {
    final Iterator<OrientedVariant> mIncluded;
    final Iterator<Variant> mExcluded;
    OrientedVariant mCurrentIncluded;
    Variant mCurrentExcluded;
    CallIterator(List<OrientedVariant> included, List<Variant> excluded) {
      mIncluded  = included.iterator();
      mExcluded = excluded.iterator();
      mCurrentIncluded = mIncluded.hasNext() ? mIncluded.next() : null;
      mCurrentExcluded = mExcluded.hasNext() ? mExcluded.next() : null;
    }

    @Override
    public boolean hasNext() {
      return mCurrentIncluded != null || mCurrentExcluded != null;
    }

    @Override
    public VariantSummary next() {
      VariantSummary result;
      if (mCurrentIncluded == null || mCurrentExcluded != null && mCurrentExcluded.getStart() < mCurrentIncluded.getStart()) {
        result = new VariantSummary(mCurrentExcluded, false, false);
        mCurrentExcluded = mExcluded.hasNext() ? mExcluded.next() : null;
      } else {
        result = new VariantSummary(mCurrentIncluded, true, mCurrentIncluded.isAlleleA());
        mCurrentIncluded = mIncluded.hasNext() ? mIncluded.next() : null;
      }

      return result;
    }
    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }

  private static class VariantSummary {
    final Variant mVariant;
    final boolean mIncluded;
    final boolean mPhase;
    VariantSummary(Variant v, boolean include, boolean alternate) {
      mVariant = v;
      mIncluded = include;
      mPhase = alternate;

    }
    boolean isPhased() {
      return mVariant.isPhased();
    }
    int startPos() {
      return mVariant.getStart();
    }
    boolean included() {
      return mIncluded;
    }
    boolean phase() {
      if (isPhased()) {
        return mPhase;
      } else {
        throw new IllegalArgumentException();
      }
    }
    @Override
    public String toString() {
      return startPos() + ":" + isPhased() + " " + (isPhased() ? (phase() ? "+" : "-") : "");
    }

  }

  @TestClass("com.rtg.variant.eval.VcfEvalTaskTest")
  static final class VariantPositionComparator implements Comparator<Variant>, Serializable {
    @Override
    public int compare(Variant o1, Variant o2) {

      if (o1.getStart() < o2.getStart()) {
        return -1;
      } else if (o1.getStart() > o2.getStart()) {
        return 1;
      }
      if (o1.getEnd() < o2.getEnd()) {
        return -1;
      } else if (o1.getEnd() > o2.getEnd()) {
        return 1;
      }
      return 0;
    }
  }

  private void merge(List<Variant> falsePositives, List<? extends Variant> b) {
    if (b.size() == 0) {
      return;
    }
    falsePositives.addAll(b);
    Collections.sort(falsePositives, new VariantPositionComparator());

  }
}
