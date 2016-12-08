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

import java.util.LinkedList;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;

/**
 * Check and mark consecutive complex calls which are equivalent.
 */
public class EquivalentFilter {

  private final byte[] mTemplateBytes;

  private List<Variant> mLastCalls;

  EquivalentFilter(byte[] refNt, List<Variant> lastCalls) {
    mTemplateBytes = refNt;
    mLastCalls = lastCalls;
  }

  /**
   * @param calls for this chunk (after complex calling but before filtering).
   * @param maxReadLength maximum read length seen in this chunk.
   * @return filtered calls (not including the last entry which needs to be held over for equivalence checking).
   */
  List<Variant> filter(List<Variant> calls, int maxReadLength) {
    // We need to be careful to check all appropriate pairs of complex calls.
    // Essentially we want to compare consecutive complex calls, ignoring most other
    // calls occurring in between (e.g. those produced by the "ALL" output mode).
    // However, interesting SNPs between complex calls sometimes break the symmetry
    // that would otherwise result in equivalent complex calls.
    // See testComplexEquivalentWithBrokenSymmetry for an example.
    final List<Variant> res = new LinkedList<>();
    for (final Variant v : calls) {
      if (!(v.isComplexScored() || v.isInteresting())) {
        if (mLastCalls == null) {
          res.add(v);
        } else {
          mLastCalls.add(v);
        }
      } else {
        if (mLastCalls != null) {
          final Variant lastCall = mLastCalls.get(0);
          if (complexEquivalent(lastCall, v, maxReadLength)) {
            v.setComplexEquivalent();
            v.addFilter(VariantFilter.COMPLEX_EQUIVALENT);
            lastCall.setComplexEquivalent();
          }
          res.addAll(mLastCalls);
        }
        mLastCalls = new LinkedList<>();
        mLastCalls.add(v);
      }
    }
    return res;
  }

  List<Variant> lastCall() {
    return mLastCalls;
  }

  private boolean complexEquivalent(Variant lastComplex, Variant current, int maxReadLength) {
    if (lastComplex == null || current == null) {
      return false;
    }
    // Efficiency, only complex calls can be equivalent
    if (!lastComplex.isComplexScored() || !current.isComplexScored()) {
      return false;
    }
    // Tests show the following does happen --- but why?
    if (lastComplex.getNumberOfSamples() != current.getNumberOfSamples()) {
      return false;
    }
    if (Math.abs(current.getLocus().getStart() - lastComplex.getLocus().getStart()) > maxReadLength) {
      return false;
    }
    boolean nonRef = false;
    for (int i = 0; i < lastComplex.getNumberOfSamples(); ++i) {
      final VariantSample ls = lastComplex.getSample(i);
      final VariantSample cs = current.getSample(i);
      if ((ls == null) != (cs == null)) { // This can happen for male model on Y chromosome spanning PAR region boundary. Not equivalent
        return false;
      }
      if (ls == null || cs == null) { // This happens regularly for female model on Y chromosome. Ignore current sample here.
        continue;
      }
      if (ls.getName() == null || cs.getName() == null) { // Overcoverage indicator pseudo-variants
        return false;
      }
      if (!ls.isIdentity() || !cs.isIdentity()) {
        nonRef = true;
      }
      if (!complexEquivalentIndividual(ls, cs, lastComplex.getLocus(), current.getLocus())) {
        return false;
      }
    }
    return nonRef; // Only variant sites can be complex equivalent
  }

  private boolean complexEquivalentIndividual(VariantSample lastComplex, VariantSample current, VariantLocus lastLocus, VariantLocus currentLocus) {
    final int startPos = lastLocus.getStart();
    final int endPos = currentLocus.getEnd();
    final String template = DnaUtils.bytesToSequenceIncCG(mTemplateBytes, startPos, endPos - startPos);
    final String[] lastPair = VariantUtils.normalizePair(lastComplex.getName());
    final String[] currentPair = VariantUtils.normalizePair(current.getName());
    return complexEquivalent(lastPair, lastLocus.getEnd() - lastLocus.getStart(), currentPair, currentLocus.getEnd() - currentLocus.getStart(), template);
  }

  /**
   * Play both calls into the template and see if they are equivalent.
   * @param lastPair the first complex call
   * @param lastLength the length of the first complex call
   * @param currentPair the second complex call
   * @param currentLength the length of the second complex call
   * @param template the template sequence from start of the first call to end of the second call
   * @return true if the calls can generate identical haplotypes
   */
  private boolean complexEquivalent(String[] lastPair, int lastLength, String[] currentPair, int currentLength, String template) {
    if (lastLength < 0) {
        return false;
    }
    final String lastTemplate = template.substring(lastLength);
    final String currentTemplate = template.substring(0, template.length() - currentLength);
    if ((lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    if ((lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    return false;
  }

}
