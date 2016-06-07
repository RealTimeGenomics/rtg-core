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

import com.rtg.variant.bayes.EvidenceAcceptor;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Model;

/**
 * This class essentially just counts the number of indels seen at a particular point, in order
 * to contribute toward triggering the complex caller.
 *
 */
public class IndelDetector implements EvidenceAcceptor {

  /** The minimum number of indels to occur at a position before invoking the complex caller */
  public static final int TRIVIAL_INDEL_COUNT = 2; //Integer.parseInt(System.getProperty("rtg.modelindel.trivial-indels", "2"));

  private int mNonTrivialInsertCount;
  private int mNonTrivialDeletionCount;
  private int mSoftClipLeftCount;
  private int mSoftClipRightCount;
  private int mNonIndelCount;

  private int mIndelLength;
  private int mSoftClipLength;

  /**
   * Constructor
   */
  public IndelDetector() {
  }

  @Override
  public void increment(EvidenceInterface evidence) {
    if (evidence == null) {
      mNonIndelCount++;
    } else {
      assert evidence instanceof EvidenceIndel;
      if (evidence.mapError() < Model.AMBIGUITY_THRESHOLD) {
        switch (evidence.read()) {
          case EvidenceIndel.INSERT:
            mNonTrivialInsertCount++;
            mIndelLength = Math.max(((EvidenceIndel) evidence).maxOperationLength(), mIndelLength);
            break;
          case EvidenceIndel.DELETE:
            mNonTrivialDeletionCount++;
            mIndelLength = Math.max(((EvidenceIndel) evidence).maxOperationLength(), mIndelLength);
            break;
          case EvidenceIndel.SOFT_CLIP_LEFT:
            mSoftClipLeftCount++;
            mSoftClipLength = Math.max(((EvidenceIndel) evidence).maxOperationLength(), mSoftClipLength);
            break;
          case EvidenceIndel.SOFT_CLIP_RIGHT:
            mSoftClipRightCount++;
            mSoftClipLength = Math.max(((EvidenceIndel) evidence).maxOperationLength(), mSoftClipLength);
            break;
          default:
            throw new RuntimeException("Invalid indel code: " + evidence.read());
        }
      }
    }
  }

  @Override
  public void freeze() {
    // Nothing to do here
  }

  /**
   * @return maximum length of an insert or delete as specified by cigar
   */
  public int maxIndelLength() {
    return mIndelLength;
  }

  /**
   *
   * @return maximum length of a soft clip as specified by cigar
   */
  public int maxSoftClipLength() {
    return mSoftClipLength;
  }

  /**
   * Get non-trivial insert count.
   * @return Returns the non-trivial insert count.
   */
  int nonTrivialInsertCount() {
    return mNonTrivialInsertCount;
  }

  /**
   * Get non-trivial deletion count.
   * @return Returns the non-trivial deletion count.
   */
  int nonTrivialDeletionCount() {
    return mNonTrivialDeletionCount;
  }

  int softClipLeftCount() {
    return mSoftClipLeftCount;
  }

  int softClipRightCount() {
    return mSoftClipRightCount;
  }

  int minIndelCount(int percentage) {
    return Math.max(TRIVIAL_INDEL_COUNT, mNonIndelCount * percentage / 100);
  }

}
