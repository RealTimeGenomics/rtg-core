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


import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;

/**
 * Maintains counts for indels from CIGAR format input.
 */
public class IndelMatcher extends EvidenceMatcher<IndelDetector> {

  /**
   * Construct a new indel matcher
   * @param template nucleotides for current sequence in template.
   * @param start of region on template to be processed (0 based, inclusive)
   */
  public IndelMatcher(final byte[] template, int start) {
    super(new ReferenceBasedBuffer<>(IndelDetectorFactory.SINGLETON, template, start), EvidenceIndelFactory.SINGLETON);
  }

  /**
   * Output matching categories and move one step in buffer.
   * @param refName name of reference sequence
   * @param startPos position in reference sequence (zero based).
   * @param endPos end position in reference sequence (zero based).
   * @param params if true then write supporting evidence and all calls even if
   * same or below threshold.
   * @return Variant object, or null if the call was not interesting.
   */
  public Variant output(String refName, int startPos, int endPos, VariantParams params) {
    final IndelDetector indelDetector = step(startPos);

    if (indelDetector != null) {
      final int minIndelCount = indelDetector.minIndelCount(params.indelTriggerFraction());
      final Ploidy ploidy = params.ploidy() == DefaultFallback.HAPLOID ? Ploidy.HAPLOID : Ploidy.DIPLOID;
      //System.err.println("@" + startPos + " nic=" + model.nonIndelCount() + " mic=" + minIndelCount + " nti=" + model.nonTrivialInsertCount() + " ntd=" + model.nonTrivialDeletionCount() + " itf=" + params.indelTriggerFraction());
      int newEnd = endPos;
      if (indelDetector.nonTrivialDeletionCount() >= minIndelCount) {
        newEnd++;
      }
      final VariantLocus locus = new VariantLocus(refName, startPos, newEnd);
      final Variant call = new Variant(locus, new VariantSample(ploidy)); //ugh... no idea what the variant sample should be in this case
      if (indelDetector.nonTrivialInsertCount() >= minIndelCount || indelDetector.nonTrivialDeletionCount() >= minIndelCount) {
        call.setInteresting();
      }
      call.setIndel(indelDetector.maxIndelLength());

      if (call.isInteresting()) {
        return call;
      }
      return null;
    } else { //reference nt == N
      return null;
    }
  }
}
