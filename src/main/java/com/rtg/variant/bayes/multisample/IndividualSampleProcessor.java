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


import com.rtg.variant.SamToMatch;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.EvidenceMatcher;
import com.rtg.variant.bayes.snp.IndelDetector;
import com.rtg.variant.bayes.snp.IndelMatcher;

/**
 * Updates underlying per-position models for every position covered by an alignment record.
 * @param <D> type of hypotheses description that this processor will handle
 */
public class IndividualSampleProcessor<D extends Description> {

  private final byte[] mTemplateBytes;
  private final EvidenceMatcher<ModelInterface<D>> mMatcherCurrent;
  private final IndelMatcher mMatcherIndel;
  private final SamToMatch mToMatch;

  /**
   * @param templateBytes nucleotides in template.
   * @param matcherCurrent matcher for incrementing models.
   * @param indelMatcher matcher for indel detection
   * @param toMatch to get models once updated.
   */
  IndividualSampleProcessor(byte[] templateBytes, EvidenceMatcher<ModelInterface<D>> matcherCurrent, IndelMatcher indelMatcher, SamToMatch toMatch) {
    mTemplateBytes = templateBytes;
    mMatcherCurrent = matcherCurrent;
    mToMatch = toMatch;
    mMatcherIndel = indelMatcher;
  }

  boolean processAlignmentRecord(final VariantAlignmentRecord var) {
    final int start = mToMatch.start(var);
    if (start < 0 || start >= mTemplateBytes.length) {
      return false;
    }
    return mToMatch.process(mTemplateBytes, var);
  }

  ModelInterface<D> step(final int start) {
    return mMatcherCurrent.step(start);
  }

  IndelDetector stepIndel(final int start) {
    return mMatcherIndel.step(start);
  }

  Variant indelOutput(String refName, int start, VariantParams params) {
    return mMatcherIndel.output(refName, start, params);
  }
}
