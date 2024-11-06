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
package com.rtg.variant.match;

import com.rtg.mode.DnaUtils;
import com.rtg.variant.NoQualityException;
import com.rtg.mode.DNARangeNAT;

/**
 * Match object where reference is returned via read methods. Nothing is
 * returned via reference methods.
 */
public class AlleleAsReadMatch extends Match {

  private final byte[] mReferenceBytes;
  private final double mQuality;
  private final int mAlleleCount;

  /**
   * Constructor
   * @param referencesBytes bytes of reference as from SDF
   * @param quality phred quality for default quality
   */
  public AlleleAsReadMatch(final byte[] referencesBytes, final double quality) {
    this(referencesBytes, quality, 0);
  }

  /**
   * Constructor
   * @param referencesBytes bytes of reference as from SDF
   */
  public AlleleAsReadMatch(final byte[] referencesBytes) {
    this (referencesBytes, 0.0, 0);
  }

  /**
   * Constructor
   * @param referencesBytes bytes of reference as from SDF
   * @param quality phred quality for default quality
   * @param alleleCount count of alleles that contributed to this match
   */
  AlleleAsReadMatch(final byte[] referencesBytes, double quality, int alleleCount) {
    mReferenceBytes = referencesBytes.clone();
    mQuality = quality;
    mAlleleCount = alleleCount;
  }

  @Override
  public boolean isFixedLeft() {
    return true;
  }

  @Override
  public boolean isFixedRight() {
    return true;
  }

  @Override
  public int length() {
    return mReferenceBytes.length;
  }

  @Override
  public double mapError() {
    return 0;
  }

  @Override
  public double baseError(final int index) throws NoQualityException {
    return mQuality;
  }

  @Override
  public double baseError() {
    return mQuality;
  }

  @Override
  public int read(final int index) {
    final int res = mReferenceBytes[index] - 1;
    assert DNARangeNAT.DNA.valid(res) : res;
    return res;
  }

  @Override
  public String readString() {
    return DnaUtils.bytesToSequenceIncCG(mReferenceBytes);
  }

  @Override
  public int alleleCount() {
    return mAlleleCount;
  }


}
