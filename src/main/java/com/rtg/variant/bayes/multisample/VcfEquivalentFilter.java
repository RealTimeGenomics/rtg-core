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

import java.util.LinkedList;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Check and mark consecutive complex calls which are equivalent.
 */
public class VcfEquivalentFilter {

  private final byte[] mTemplateBytes;

  private List<VcfRecord> mLastCalls;

  VcfEquivalentFilter(byte[] refNt, List<VcfRecord> lastCalls) {
    mTemplateBytes = refNt;
    mLastCalls = lastCalls;
  }

  static boolean isComplexScored(final VcfRecord rec) {
    return rec.hasInfo("XRX");
  }

  /**
   * @param calls for this chunk (after complex calling but before filtering).
   * @param maxReadLength maximum read length seen in this chunk.
   * @return filtered calls (not including the last entry which needs to be held over for equivalence checking).
   */
  List<VcfRecord> filter(List<VcfRecord> calls, int maxReadLength) {
    // We need to be careful to check all appropriate pairs of complex calls.
    // Essentially we want to compare consecutive complex calls, ignoring most other
    // calls occurring in between (e.g. those produced by the "ALL" output mode).
    // However, interesting SNPs between complex calls sometimes break the symmetry
    // that would otherwise result in equivalent complex calls.
    // See testComplexEquivalentWithBrokenSymmetry for an example.
    final List<VcfRecord> res = new LinkedList<>();
    for (final VcfRecord v : calls) {
      if (!isComplexScored(v)) {
        if (mLastCalls == null) {
          res.add(v);
        } else {
          mLastCalls.add(v);
        }
      } else {
        if (mLastCalls != null && !mLastCalls.isEmpty()) {
          final VcfRecord lastCall = mLastCalls.get(0);
          if (complexEquivalent(lastCall, v, maxReadLength)) {
            v.addFilter("RCEQUIV");
            v.setInfo("RCE");
            lastCall.setInfo("RCE");
          }
          res.addAll(mLastCalls);
        }
        mLastCalls = new LinkedList<>();
        mLastCalls.add(v);
      }
    }
    return res;
  }

  List<VcfRecord> lastCall() {
    return mLastCalls;
  }

  private boolean complexEquivalent(VcfRecord lastComplex, VcfRecord current, int maxReadLength) {
    if (lastComplex == null || current == null) {
      return false;
    }
    // Efficiency, only complex calls can be equivalent
    if (!isComplexScored(lastComplex) || !isComplexScored(current)) {
      return false;
    }
    // Tests show the following does happen --- but why?
    if (lastComplex.getNumberOfSamples() != current.getNumberOfSamples()) {
      return false;
    }
    if (!current.getSequenceName().equals(lastComplex.getSequenceName())) {
      return false;
    }
    if (Math.abs(current.getStart() - lastComplex.getStart()) > maxReadLength) {
      return false;
    }
    boolean nonRef = false;
    for (int i = 0; i < lastComplex.getNumberOfSamples(); ++i) {
      final int[] ls = VcfUtils.getValidGt(lastComplex, i);
      final int[] cs = VcfUtils.getValidGt(current, i);
      if (VcfUtils.isMissingGt(ls) != VcfUtils.isMissingGt(cs)) {
        return false; // This can happen for male model on Y chromosome spanning PAR region boundary. Not equivalent
      }
      if (VcfUtils.isMissingGt(ls) || VcfUtils.isMissingGt(cs)) {
        continue; // This happens regularly for female model on Y chromosome. Ignore current sample here.
      }
      if (ls.length != cs.length) {
        return false; // Can't be equivalent if ploidy is different
      }
      if (!VcfUtils.isHomozygousRef(ls) || !VcfUtils.isHomozygousRef(cs)) {
        nonRef = true;
      }
      final int startPos = lastComplex.getStart();
      final int endPos = current.getEnd();
      final String template = DnaUtils.bytesToSequenceIncCG(mTemplateBytes, startPos, endPos - startPos);
      final String[] lastPair = getAlleles(lastComplex, ls);
      final String[] currentPair = getAlleles(current, cs);
      if (!complexEquivalent(lastPair, lastComplex.getEnd() - startPos, currentPair, endPos - current.getStart(), template)) {
        return false;
      }
    }
    return nonRef; // Only variant sites can be complex equivalent
  }

  private static String[] getAlleles(final VcfRecord rec, final int[] gt) {
    final String[] ret = new String[gt.length];
    int i = 0;
    for (final int allele : gt) {
      ret[i++] = allele == 0 ? rec.getRefCall() : rec.getAltCalls().get(allele - 1);
    }
    return ret;
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
    if (lastLength < 0 || lastPair.length != currentPair.length || lastPair.length < 1) {
        return false;
    }
    final String lastTemplate = template.substring(lastLength);
    final String currentTemplate = template.substring(0, template.length() - currentLength);
    if (lastPair.length == 1) {
      return (lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[0]);
    }
    // Check both phasings - note this won't work for higher ploidys
    if ((lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    if ((lastPair[1] + lastTemplate).equals(currentTemplate + currentPair[0]) && (lastPair[0] + lastTemplate).equals(currentTemplate + currentPair[1])) {
      return true;
    }
    return false;
  }

}
