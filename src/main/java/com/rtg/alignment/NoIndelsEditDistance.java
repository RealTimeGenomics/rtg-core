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
package com.rtg.alignment;

import com.rtg.mode.DNA;
import com.rtg.ngs.NgsParams;

/**
 * Edit distance which only aligns mismatches, allowing as many as possible before an
 * indel could be preferred.
 * Unknown nucleotides are unhandled.
 */
class NoIndelsEditDistance implements UnidirectionalEditDistance {

  private static final int UNKNOWN = DNA.N.ordinal();

  private int[] mWorkspace;

  private final int[] mMismatchPositions;

  private final int mSubstitutionPenalty;

  private int mMaxReadLengthSeen = 0;

  NoIndelsEditDistance(NgsParams ngsParams) {

    /*
     * We want to allow this aligner to account for all alignments which can be expressed with no indels.
     * A single indel gives an alignment penalty of (gapOpen + gapExtend), so we want as many substitutions
     * which can "fit" into that alignment score. We prefer substitutions to indels, so
     * if a score exists for example:
     * if gapOpen + gapExtend = 3, and substitutionPenalty = 1
     * We would prefer 3 substitutions than one indel, and this aligner should accept 3 or fewer mismatches.
     * However, if the substitutionPenalty is 2, then we can only safely find 1 mismatch.
     */

    final int maxSubstitutions = (ngsParams.gapOpenPenalty() + ngsParams.gapExtendPenalty()) / ngsParams.substitutionPenalty();
    mMismatchPositions = new int[maxSubstitutions];
    mSubstitutionPenalty = ngsParams.substitutionPenalty();
  }

  @Override
  public void logStats() {
    // do nothing
  }

  /**
   * Tries to perform a quick alignment with substitutions only and as many
   * as possible, permitted by the penalties. This is not
   * useful for Complete Genomics data (and all Complete Genomics testing has
   * been removed from it).
   *
   * @param rlen length of read (assumed that given read is at least this long)
   * @param template the encoded template
   * @param zeroBasedStart start position in the template
   * @param maxScore currently ignored.
   * @param maxShift ignored. No shifting permitted.
   * @param cgLeft ignored.
   * @return the actions array, or null if not possible
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    // Quick idiocy checks first to avoid string construction cost
    if (zeroBasedStart < 0 || zeroBasedStart + rlen > template.length) {    //if read mapped off template at either end, just delegate to next in chain
      return null;
    }
    int mismatchIndex = 0; // Mismatch count, fail if it exceeds mMismatchPositions.length
    // Our use of multiplication to quickly test for UNKNOWNs in the following
    // code requires UNKNOWN to always equal zero (as documented in DNA.java).
    assert UNKNOWN == 0;

    //ensure workspace is big enough
    if (rlen > mMaxReadLengthSeen) {
      final int size = ActionsHelper.ACTIONS_START_INDEX + 1 + (int) (rlen / (double) ActionsHelper.ACTIONS_PER_INT + 0.5);
      if (mWorkspace == null || mWorkspace.length < size) {
        mWorkspace = new int[size];
      }
      mMaxReadLengthSeen = rlen;
    }

    for (int readPos = 0, tPos = zeroBasedStart; readPos < rlen && tPos < template.length; ++readPos, ++tPos) {
      final byte t = template[tPos];
      final byte r = read[readPos];
      // Note: the multiplication ignores mismatches when one/both are UNKNOWN.

      if (r * t == UNKNOWN) {
        return null;
      }
      if (r != t) {   //nt different
        if (mismatchIndex >= mMismatchPositions.length) {
          return null;    //we've exceeded the substitution threshold, delegate to next in chain
        } else if ((mismatchIndex + 1) * mSubstitutionPenalty > maxScore) { //TODO test this stuff by overriding the earlyTerminate method in tests.
          return earlyTerminate(zeroBasedStart);
        }

        mMismatchPositions[mismatchIndex] = readPos;
        mismatchIndex += 1;
      }
    }

    // Okay, can build a quick solution
    mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = 0;
    mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;

    if (mismatchIndex == 0) {
      ActionsHelper.prepend(mWorkspace, rlen, ActionsHelper.SAME, 0); //no mismatches so just return perfect array
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
      mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = rlen;
    } else {
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart + rlen; //add rlen because the actionshelper prepend "helpfully" alters the workspace start position as it goes
      int previousMismatchPos = rlen;  //initially say there was a mismatch at rlen
      for (int i = mismatchIndex; i > 0; --i) { //if mismatchIndex is 0, then we haven't found any mismatches, so don't evaluate that.
        final int thisMismatchPos = mMismatchPositions[i - 1];
        final int matchLengthToOutput = previousMismatchPos - thisMismatchPos - 1;

        if (matchLengthToOutput > 0) {
          ActionsHelper.prepend(mWorkspace, matchLengthToOutput, ActionsHelper.SAME, 0); //no mismatches so just return perfect array
        }
        ActionsHelper.prepend(mWorkspace, 1, ActionsHelper.MISMATCH, mSubstitutionPenalty);
        previousMismatchPos = thisMismatchPos;
      }

      //output final chunk of matches if necessary
      if (previousMismatchPos > 0) {
        ActionsHelper.prepend(mWorkspace, previousMismatchPos, ActionsHelper.SAME, 0);
      }
    }
    return mWorkspace;
  }

  private int[] earlyTerminate(int zeroBasedStart) {
    mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
    mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
    mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
    return mWorkspace;
  }
  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos,
      int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException(); //no good reason this couldn't be supported if desired (as long as length of read and length of template are same)
  }
  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos,
      int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException(); //no good reason this couldn't be supported if desired.
  }
  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException(); //no good reason this couldn't be supported if desired.
  }
}
