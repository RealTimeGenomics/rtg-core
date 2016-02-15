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
package com.rtg.variant.bayes.snp;


import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;

/**
 * Maintains counts for indels from CIGAR format input.
 */
public class IndelMatcher extends EvidenceMatcher<IndelDetector> {

  /**
   * Construct a new indel matcher
   * @param buffer manages per-position <code>EvidenceAcceptor</code> s
   */
  public IndelMatcher(ReferenceBasedBuffer<IndelDetector> buffer) {
    super(buffer, EvidenceIndelFactory.SINGLETON);
  }

  /**
   * Output matching categories and move one step in buffer.
   * @param refName name of reference sequence
   * @param startPos position in reference sequence (zero based).
   * @param params if true then write supporting evidence and all calls even if
   * same or below threshold.
   * @return Variant object, or null if the call was not interesting.
   */
  public Variant output(String refName, int startPos, VariantParams params) {
    final IndelDetector indelDetector = step(startPos);

    if (indelDetector != null) {
      final int minIndelCount = indelDetector.minIndelCount(params.indelTriggerFraction());
      //System.err.println("@" + startPos + " nic=" + model.nonIndelCount() + " mic=" + minIndelCount + " nti=" + model.nonTrivialInsertCount() + " ntd=" + model.nonTrivialDeletionCount() + " itf=" + params.indelTriggerFraction());
      int newEnd = startPos;
      if (indelDetector.nonTrivialDeletionCount() >= minIndelCount) {
        newEnd++;
      }
      if (indelDetector.nonTrivialInsertCount() >= minIndelCount || indelDetector.nonTrivialDeletionCount() >= minIndelCount) {
        final VariantLocus locus = new VariantLocus(refName, startPos, newEnd);
        final Variant call = new Variant(locus);
        call.setInteresting();
        call.setIndel(indelDetector.maxIndelLength());
        return call;
      }
      return null;
    } else { //reference nt == N
      return null;
    }
  }
}
