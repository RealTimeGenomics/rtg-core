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
    return mToMatch.process(mTemplateBytes, var, null);
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
