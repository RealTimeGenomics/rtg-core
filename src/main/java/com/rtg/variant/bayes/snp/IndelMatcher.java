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
package com.rtg.variant.bayes.snp;


import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;

/**
 * Maintains counts for indels from CIGAR format input.
 */
public class IndelMatcher extends EvidenceMatcher<IndelDetector> {

  private static final boolean DUMP_INDEL_SIGNAL = GlobalFlags.getBooleanValue(CoreGlobalFlags.DUMP_COMPLEX_TRIGGER_SIGNALS);

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
        ++newEnd;
      }
      if (DUMP_INDEL_SIGNAL) {
        if (indelDetector.nonTrivialDeletionCount() > 1 || indelDetector.nonTrivialInsertCount() > 1 || indelDetector.softClipLeftCount() > 1 || indelDetector.softClipRightCount() > 1) {
          Diagnostic.developerLog(String.format("INDEL-SIGNAL\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
            refName, startPos, newEnd,
            indelDetector.totalCount(),
            indelDetector.nonTrivialInsertCount(), indelDetector.nonTrivialDeletionCount(),
            indelDetector.softClipLeftCount(), indelDetector.softClipRightCount(),
            indelDetector.indelLength(), indelDetector.maxSoftClipLength()));
        }
      }
      if (indelDetector.nonTrivialInsertCount() >= minIndelCount || indelDetector.nonTrivialDeletionCount() >= minIndelCount) {
        final VariantLocus locus = new VariantLocus(refName, startPos, newEnd);
        final Variant call = new Variant(locus);
        call.setInteresting();
        call.setIndel(indelDetector.indelLength());
        return call;
      }
      if (indelDetector.softClipLeftCount() >= minIndelCount || indelDetector.softClipRightCount() >= minIndelCount) {
        final VariantLocus locus = new VariantLocus(refName, startPos, newEnd);
        final Variant call = new Variant(locus);
        call.setInteresting();
        final Variant.SoftClipSide side = indelDetector.softClipLeftCount() > indelDetector.softClipRightCount() ? Variant.SoftClipSide.LEFT : Variant.SoftClipSide.RIGHT;
        call.setSoftClip(indelDetector.maxSoftClipLength(), side);
        return call;
      }
      return null;
    } else { //reference nt == N
      return null;
    }
  }
}
