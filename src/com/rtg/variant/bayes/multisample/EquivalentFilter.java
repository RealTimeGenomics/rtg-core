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
//    System.err.println("last: " + lastComplex);
//    System.err.println("current: " + current);
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

    //We might have to consider top level call in deciding this
    for (int i = 0; i < lastComplex.getNumberOfSamples(); i++) {
      if (!complexEquivalentIndividual(lastComplex.getSample(i), current.getSample(i), lastComplex.getLocus(), current.getLocus(), maxReadLength)) {
        return false;
      }
    }
    return true;
  }

  private boolean complexEquivalentIndividual(VariantSample lastComplex, VariantSample current, VariantLocus lastLocus, VariantLocus currentLocus, int maxReadLength) {
    //System.err.println("lastComplex=" + lastComplex);
    //System.err.println("current=" + current);
    if (lastComplex == null || current == null) {
      return false;
    }
//    if (!lastComplex.isAllowOutput() || !current.isAllowOutput()) {
//      return false;
//    }
    if (Math.abs(currentLocus.getStart() - lastLocus.getStart()) > maxReadLength) {
      return false;
    }
    if (lastComplex.getName() == null || current.getName() == null) {
      return false;
    }
    return checkComplexEquivalent(lastComplex, current, lastLocus, currentLocus);
  }

  private boolean checkComplexEquivalent(VariantSample lastComplex, VariantSample current, VariantLocus lastLocus, VariantLocus currentLocus) {
    final int startPos = lastLocus.getStart();
    final int endPos = currentLocus.getEnd();
    final String template = DnaUtils.bytesToSequenceIncCG(mTemplateBytes, startPos, endPos - startPos);
    final String[] lastPair = VariantUtils.normalizePair(lastComplex.getName());
    final String[] currentPair = VariantUtils.normalizePair(current.getName());
    // System.err.println("**" + lastComplex);
    // System.err.println("name=" + lastComplex.getName() + " currentName=" + current.getName());
    // System.err.println("startPos=" + startPos +" endPos=" + endPos);
    // Using getOutputPosition() here makes this safe to do post-trimming
    return complexEquivalent(lastPair, lastLocus.getEnd() - lastLocus.getStart(), currentPair, currentLocus.getEnd() - currentLocus.getStart(), template);
  }

  /**
   * @param lastPair the first complex call
   * @param lastLength the length of the first complex call
   * @param currentPair the second complex call
   * @param currentLength the length of the second complex call
   * @param template the template sequence
   * @return true if the calls generate an identical template
   */
  private boolean complexEquivalent(String[] lastPair, int lastLength, String[] currentPair, int currentLength, String template) {
    if (lastLength < 0) {
        return false;
    }
    //System.err.println("lastLength=" + lastLength + " template=\"" + template + "\"");
    final String lastTemplate = template.substring(lastLength);
    final String currentTemplate = template.substring(0, template.length() - currentLength);
    //    System.err.println(lastTemplate + " " + currentTemplate);
    //    System.err.println(lastPair[0] + lastTemplate + ":" + currentTemplate + currentPair[0]);
    //    System.err.println(lastPair[1] + lastTemplate + ":" + currentTemplate + currentPair[1]);
    if ((lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    if ((lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    return false;
  }

}
