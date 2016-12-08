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

/**
 */
public class SoftClipperOmni implements EditDistance {

  private final int mClipIndelsXFromEnd;

  private final EditDistance mEd;

  private final ActionsHelper.CommandIterator mForwardIterator;
  private final ActionsHelper.CommandIterator mReverseIterator;

  private int mStartPositionOffset;
  private int mDels;

  /**
   * Create an edit distance wrapper which soft clips indels near the ends of alignments.
   * @param ed an edit distance implementation to wrap
   * @param clipIndelsXFromEnd the number of bases from each end of an alignment to check for indels
   */
  public SoftClipperOmni(EditDistance ed, int clipIndelsXFromEnd) {
    mEd = ed;
    mClipIndelsXFromEnd = clipIndelsXFromEnd;
    mForwardIterator = new ActionsHelper.CommandIterator();
    mReverseIterator = new ActionsHelper.CommandIteratorReverse();
  }

  int[] softClipActions(int[] actions, boolean rc) {
    mReverseIterator.setActions(actions);
    int softClips = numSoftClips(mReverseIterator);
    if (softClips > 0) {
      ActionsHelper.softClip(actions, false, softClips - mDels, mDels);
      if (rc) {
        ActionsHelper.setZeroBasedTemplateStart(actions, ActionsHelper.zeroBasedTemplateStart(actions) + softClips + mStartPositionOffset);
      }
    }

    mForwardIterator.setActions(actions);
    softClips = numSoftClips(mForwardIterator);
    if (softClips > 0) {
      ActionsHelper.softClip(actions, true, softClips - mDels, mDels);
      if (!rc) {
        ActionsHelper.setZeroBasedTemplateStart(actions, ActionsHelper.zeroBasedTemplateStart(actions) + softClips + mStartPositionOffset);
      }
    }

    return actions;
  }

  /**
   * Detects INDELs within <code>mClipIndelsXFromEnd</code> of the "start" of an actions command iterator.
   * Also sets <code>mStartPositionOffset</code> to a value which can be added to the start position to account for this.
   * @param it a <code>CommandIterator</code> to iterate across.
   * @return the number of actions which are to be soft clipped.
   */
  private int numSoftClips(ActionsHelper.CommandIterator it) {
    int pos = 0;
    int softClipPos = -1;
    mStartPositionOffset = 0;
    mDels = 0;
    boolean extend = false;
    while (it.hasNext() && pos < mClipIndelsXFromEnd) {
      final int command = it.next();
      switch (command) {
        case ActionsHelper.INSERTION_INTO_REFERENCE:
          softClipPos = pos;
          --mStartPositionOffset;
          extend = true;
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          softClipPos = pos;
          ++mDels;
          extend = true;
          break;
        default:
          if (extend && command != ActionsHelper.SAME) {
            ++softClipPos;
          } else {
            extend = false;
          }
      }
      ++pos;
    }

    //if there's an indel continuing past the cutoff position, we'll need to soft clip that too.
    if (extend) {
      while (it.hasNext()) {
        final int command = it.next();
        if (command != ActionsHelper.SAME) {
          switch (command) {
            case ActionsHelper.INSERTION_INTO_REFERENCE:
              --mStartPositionOffset;
              break;
            case ActionsHelper.DELETION_FROM_REFERENCE:
              ++mDels;
              break;
            default:
          }
          ++softClipPos;
        } else {
          break;
        }
      }
    }
    return softClipPos == -1 ? 0 : softClipPos + 1;
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
    final int[] res = mEd.calculateEditDistance(read, rlen, template, zeroBasedStart, rc, maxScore, maxShift, cgLeft);
    if (mClipIndelsXFromEnd <= 0 || res == null || res[ActionsHelper.ALIGNMENT_SCORE_INDEX] == Integer.MAX_VALUE) {
      return res;
    }
    return softClipActions(res, rc);
  }

  @Override
  public void logStats() {
    mEd.logStats();
  }
}
