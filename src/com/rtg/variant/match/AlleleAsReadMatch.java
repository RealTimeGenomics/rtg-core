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
package com.rtg.variant.match;

import com.rtg.mode.DnaUtils;
import com.rtg.variant.NoQualityException;
import com.rtg.variant.dna.DNARangeNAT;

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
