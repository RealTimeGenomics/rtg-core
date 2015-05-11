/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 */
public class PhasingEvaluator {

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


  static PhasingResult countMisphasings(Path best) {
    final CallIterator baseline = new CallIterator(best.getBaselineIncluded(), best.getBaselineExcluded());
    final CallIterator calls = new CallIterator(best.getCalledIncluded(), best.getCalledExcluded());
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

      // Collect all baseline variants in the sync region
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

      // Collect all called variants in the sync region
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
      final VariantSummary result;
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
}
