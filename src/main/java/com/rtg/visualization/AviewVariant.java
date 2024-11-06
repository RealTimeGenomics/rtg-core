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
package com.rtg.visualization;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;


/**
 * Class for holding information parsed from a SNP detection file
 */
public class AviewVariant {

  private final int mReferenceLength;
  private final int mPosition;
  private final boolean mIsFiltered;
  private final String[] mPrediction;

  /**
   * Construct a Detected Variant by inspecting a <code>VcfRecord</code> object. A string is supplied as the output representation of this object
   * @param rec VCF record to process
   * @param sampleNo the sample column number (starting from 0) for multiple sample variant calls
   */
  public AviewVariant(VcfRecord rec, int sampleNo) {
    final String gt = rec.getFormat(VcfUtils.FORMAT_GENOTYPE).get(sampleNo);
    final int[] gtArray = VcfUtils.splitGt(gt);
    if (gtArray.length > 2) {
      throw new NoTalkbackSlimException("Cannot handle calls with ploidy more than diploid.");
    }

    final boolean hasPreviousNt = VcfUtils.hasRedundantFirstNucleotide(rec);
    final String[] alleleStrings = VcfUtils.getAlleleStrings(rec, hasPreviousNt);
    mPrediction = new String[VcfUtils.isHomozygousAlt(rec, sampleNo) ? 1 : 2];
    for (int i = 0; i < mPrediction.length; ++i) {
      final int alleleId = gtArray[i];
      mPrediction[i] = alleleId == -1 ? "N" : alleleStrings[alleleId];
    }
    mIsFiltered = rec.isFiltered();
    mReferenceLength = rec.getRefCall().length() - (hasPreviousNt ? 1 : 0);
    mPosition = rec.getOneBasedStart() + (hasPreviousNt ? 1 : 0);
  }


  /**
   * @return the position (1 based).
   */
  public int getPosition() {
    return mPosition;
  }

  /**
   * @return the length of the reference section
   */
  public int referenceLength() {
    return mReferenceLength;
  }

  /**
   * One allele of the variant as determined by allele parameter.
   * @param alleleA if true select the A allele.
   * @return the allele (may be null or zero length)
   */
  public String nt(boolean alleleA) {
    if (mPrediction.length == 2) {
      return mPrediction[alleleA ? 0 : 1];
    } else if (alleleA) {
      return mPrediction[0];
    } else {
      return null;
    }
  }

  /**
   * The first allele of the variant. Always non null but may be zero length.
   * @return the A allele
   */
  public String ntAlleleA() {
    return nt(true);
  }

  /**
   * Second allele of the variant. Will be null if homozygous and may be zero length.
   * @return the B allele, or null if homozygous
   */
  public String ntAlleleB() {
    return nt(false);
  }

  /**
   * @return true if original record was filtered, false otherwise
   */
  public boolean isFiltered() {
    return mIsFiltered;
  }
}
