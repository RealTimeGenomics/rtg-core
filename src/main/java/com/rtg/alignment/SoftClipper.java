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
package com.rtg.alignment;

/**
 * Performs soft-clipping of an alignment.
 * If an indel occurs within N bases of the end of the alignment, it will be clipped.
 * If mismatches occur within M bases of the end of the alignment, it it will be clipped.
 * If the post-clipping alignment has fewer than O matches, return poorest alignment score.
 */
public class SoftClipper implements UnidirectionalEditDistance {

  private final UnidirectionalEditDistance mEd;
  private final int mMinMatches;
  private final int mClipIndelsXFromEnd;
  private final int mClipMismatchesXFromEnd;

  private final ActionsHelper.CommandIterator mForwardIterator;
  private final ActionsHelper.CommandIterator mReverseIterator;

  private int mStartPositionOffset;
  private int mDels;


  /**
   * Create an edit distance wrapper which soft clips indels and mismatches near the ends of alignments.
   * @param ed an edit distance implementation to wrap
   * @param clipIndelsXFromEnd the number of bases from each end of an alignment to check for indels
   * @param clipMismatchesXFromEnd the number of bases from each end of an alignment to check for mismatches to soft-clip
   */
  public SoftClipper(UnidirectionalEditDistance ed, int clipIndelsXFromEnd, int clipMismatchesXFromEnd) {
    this(ed, clipIndelsXFromEnd, clipMismatchesXFromEnd, 0);
  }

  /**
   * Create an edit distance wrapper which soft clips indels and mismatches near the ends of alignments.
   * @param ed an edit distance implementation to wrap
   * @param clipIndelsXFromEnd the number of bases from each end of an alignment to check for indels
   * @param clipMismatchesXFromEnd the number of bases from each end of an alignment to check for mismatches to soft-clip
   * @param minMatches the minimum number of matches required in a post-clipped alignment for the alignment to be retained
   */
  public SoftClipper(UnidirectionalEditDistance ed, int clipIndelsXFromEnd, int clipMismatchesXFromEnd, int minMatches) {
    mEd = ed;
    mClipIndelsXFromEnd = clipIndelsXFromEnd;
    mClipMismatchesXFromEnd = clipMismatchesXFromEnd;
    mMinMatches = minMatches;
    mForwardIterator = new ActionsHelper.CommandIterator();
    mReverseIterator = new ActionsHelper.CommandIteratorReverse();
  }

  int[] softClipActions(boolean start, boolean end, int[] actions) {
    if (actions == null || actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] == Integer.MAX_VALUE) {
      return actions;
    }
    if (mClipIndelsXFromEnd > 0 || mClipMismatchesXFromEnd > 0) {
      if (end) {
        mReverseIterator.setActions(actions);
        final int softClips = numSoftClips(mReverseIterator);
        if (softClips > 0) {
          ActionsHelper.softClip(actions, false, softClips - mDels, mDels);
        }
      }

      if (start) {
        mForwardIterator.setActions(actions);
        final int softClips = numSoftClips(mForwardIterator);
        if (softClips > 0) {
          ActionsHelper.softClip(actions, true, softClips - mDels, mDels);
          ActionsHelper.setZeroBasedTemplateStart(actions, ActionsHelper.zeroBasedTemplateStart(actions) + softClips + mStartPositionOffset);
        }
      }
    }
    if (mMinMatches > 0 && ActionsHelper.matchCount(actions) < mMinMatches) {
      actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
    }
    return actions;
  }

  /**
   * Determines how many bases to soft clip from the "start" of an actions command iterator.
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
    int anchor = 0;
    final int maxScan = Math.max(mClipIndelsXFromEnd, mClipMismatchesXFromEnd);
    while (pos < maxScan && it.hasNext()) {
      final int command = it.next();
      switch (command) {
        case ActionsHelper.INSERTION_INTO_REFERENCE:
          --mStartPositionOffset;
          if (extend || pos < mClipIndelsXFromEnd) {
            softClipPos = pos;
            extend = true;
          }
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          ++mDels;
          if (extend || pos < mClipIndelsXFromEnd) {
            softClipPos = pos;
            extend = true;
          }
          break;
        case ActionsHelper.SAME:
          ++anchor;
          extend = false;
          break;
        default:
          if (extend) {
            ++softClipPos;
          } else if (command == ActionsHelper.MISMATCH && anchor < mClipMismatchesXFromEnd) {
            softClipPos = pos;
            extend = true;
          }
          break;
      }
      ++pos;
    }

    // Triggered, now extend the soft clipping up to the next match if need be
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
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    return softClipActions(true, true, mEd.calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, cgLeft));
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    return softClipActions(false, true, mEd.calculateEditDistanceFixedStart(read, readStartPos, readEndPos, template, templateStartPos, maxScore, maxShift));
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
    return softClipActions(true, false, mEd.calculateEditDistanceFixedEnd(read, readStartPos, readEndPos, template, templateExpectedStartPos, templateEndPos, maxScore, maxShift));
  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
    return mEd.calculateEditDistanceFixedBoth(read, readStartPos, readEndPos, template, templateStartPos, templateEndPos, maxScore, maxShift);
  }

  @Override
  public void logStats() {
    mEd.logStats();
  }
}
