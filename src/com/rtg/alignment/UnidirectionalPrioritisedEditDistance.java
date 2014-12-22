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

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.launcher.GlobalFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.util.License;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * A prioritised edit distance which attempts multiple other edit distances in
 * order of simplicity, returning the result from the first which passes.
 *
 */
class UnidirectionalPrioritisedEditDistance implements UnidirectionalEditDistance {

  private static final boolean LOG_AS_HISTO = GlobalFlags.isSet(GlobalFlags.EDIT_DIST_LOG_AS_HISTOGRAM_FLAG);

  private final UnidirectionalEditDistance[] mEds;
  private final int[] mCounts;
  private final long[] mTimeTaken;
  private final int[] mNullReturned;
  private final int[] mMaxIntReturned;
  private final long[] mTotalScore;

  private final int[] mMaxIntActions;

  private final int[] mASHistogramUnderMaxScore;
  private final int[] mASHistogramOverMaxScore;

  private final NgsParams mNgsParams;

  private final boolean mValidateAlignments;

  /**
   * Creates a Prioritised edit distance which tries several options in order
   * @param ngsParams {@link NgsParams} for current run
   * @param editDistances the edit distances to iterate over (in order)
   */
  protected UnidirectionalPrioritisedEditDistance(NgsParams ngsParams, UnidirectionalEditDistance... editDistances) {
    mEds = editDistances;
    mCounts = new int[editDistances.length];
    mTimeTaken = new long[editDistances.length]; // scaled to be in micro seconds, estimated by sampling
    mNullReturned = new int[editDistances.length];
    mMaxIntReturned = new int[editDistances.length];
    mTotalScore = new long[editDistances.length];

    mMaxIntActions = new int[12];
    mMaxIntActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;

    mNgsParams = ngsParams;

    if (LOG_AS_HISTO) {
      mASHistogramUnderMaxScore = new int[256];  //TODO arbitrary limitation, should be based on max scores (but if those are percentages, depends on read lengths)
      mASHistogramOverMaxScore = new int[mASHistogramUnderMaxScore.length];
    } else {
      mASHistogramUnderMaxScore = null;
      mASHistogramOverMaxScore = null;
    }

    mValidateAlignments = License.isDeveloper() && "".equals(GlobalFlags.getStringValue(GlobalFlags.EDIT_DIST_INDEL_TABLE_FLAG)); //TODO: table can't currently be replicated by the validator, so don't bother checking.
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    if (rlen < 0 || zeroBasedStart >= template.length) {
      Diagnostic.developerLog("PROBLEM: UnidirectionalED problem: " + rlen + " " + zeroBasedStart + "/" + template.length + " " + maxScore + " "
          + DnaUtils.bytesToSequenceIncCG(read));
    }

    int i = 0;
    long starttime = 0;
    long endtime;
    for (final UnidirectionalEditDistance ed : mEds) {
      mCounts[i]++;
      if ((mCounts[i] & 1023) == 0) {
        starttime = System.nanoTime();
      }
      final int[] actions = ed.calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
      if ((mCounts[i] & 1023) == 0) {
        endtime = System.nanoTime();
        mTimeTaken[i] += endtime - starttime;
      }
      if (actions != null) {
        final int score = actions[ActionsHelper.ALIGNMENT_SCORE_INDEX];
        if (score == Integer.MAX_VALUE) {
          mMaxIntReturned[i]++;
          if (LOG_AS_HISTO) {
            mASHistogramUnderMaxScore[mASHistogramUnderMaxScore.length - 1]++;
          }
        } else {
         if (mValidateAlignments) {
            checkResult(read, rlen, template, zeroBasedStart, maxScore, i, actions);
            if (LOG_AS_HISTO) {
              if (score > mASHistogramOverMaxScore.length - 2) {
                assert maxScore < mASHistogramOverMaxScore.length;
                mASHistogramOverMaxScore[mASHistogramOverMaxScore.length - 2]++;
              } else {
                if (score > maxScore) {
                  mASHistogramOverMaxScore[score]++;
                } else {
                  mASHistogramUnderMaxScore[score]++;
                }
              }
            }
          }
          mTotalScore[i] += score;
        }
        return actions;
      }
      mNullReturned[i]++;
      i++;
    }
    mMaxIntActions[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
    return mMaxIntActions;
  }

  @JumbleIgnore
  private void checkResult(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int i, int[] actions) {
    // good result, check every few
    if (((mCounts[i] & ((2L << 18) - 1)) == 0) && read.length > 0 && template.length > 0) {
      // check validation
      try {
        final ActionsValidator av = new ActionsValidator(mNgsParams.gapOpenPenalty(), mNgsParams.gapExtendPenalty(), mNgsParams.substitutionPenalty(), mNgsParams.unknownsPenalty());
        if (!av.isValid(actions, read, rlen, template, maxScore)) {
          final String msg = av.getErrorDetails(actions, read, rlen, template, zeroBasedStart);
          Diagnostic.developerLog("UnidirectionalPrioritisedEditDistance validation action problem: returned by " + mEds[i].getClass().getName()
              + StringUtils.LS + msg);
        }
      } catch (final RuntimeException e) {
        Diagnostic.developerLog("UnidirectionalPrioritisedEditDistance: Problem even validating this pair of reads" + StringUtils.LS);
        Diagnostic.developerLog(" read:  " + DnaUtils.bytesToSequenceIncCG(read, 0, rlen) + StringUtils.LS);
        Diagnostic.developerLog(" tmpl:  " + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart, rlen) + StringUtils.LS);
        Diagnostic.developerLog(e);
      }
    }

  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
    for (final UnidirectionalEditDistance ed : mEds) {
      final int[] actions = ed.calculateEditDistanceFixedBoth(read, readStartPos, readEndPos, template, templateStartPos, templateEndPos, maxScore, maxShift);
      if (actions != null) {
        return actions;
      }
    }
    return mMaxIntActions;
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
    for (final UnidirectionalEditDistance ed : mEds) {
      final int[] actions = ed.calculateEditDistanceFixedEnd(read, readStartPos, readEndPos, template, templateExpectedStartPos, templateEndPos,
          maxScore, maxShift);
      if (actions != null) {
        return actions;
      }
    }
    return mMaxIntActions;
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateStartPos, int maxScore, int maxShift) {
    for (final UnidirectionalEditDistance ed : mEds) {
      final int[] actions = ed.calculateEditDistanceFixedStart(read, readStartPos, readEndPos, template, templateStartPos, maxScore, maxShift);
      if (actions != null) {
        return actions;
      }
    }
    return mMaxIntActions;
  }

  @Override
  @JumbleIgnore
  public void logStats() {
    final StringBuilder sb = new StringBuilder();
    sb.append("UniDirectionalPrioritisedEditDistance stats").append(StringUtils.LS);
    for (int i = 0; i < mCounts.length; i++) {
      final double timetaken = mTimeTaken[i] * 1024 / 1000000000.0;
      final double speed = MathUtils.round(mCounts[i] / timetaken);
      final int consumed = mCounts[i] - mNullReturned[i];
      final double eatingspeed = MathUtils.round(consumed / timetaken);
      sb.append("EDstat  [ ").append(i + 1).append(" ] calls=").append(mCounts[i]).append(" nulls=").append(mNullReturned[i]).append(" MaxInt=").append(mMaxIntReturned[i]).append(" avg.score=").append(mCounts[i] != 0 ? "" + (mTotalScore[i] / mCounts[i]) : "-").append(" time(s)=").append(timetaken).append(" speed(per s)=").append(speed).append(" eating(per s)=").append(eatingspeed).append("   ").append(mEds[i].getClass().getName());
      sb.append(StringUtils.LS);
    }
    Diagnostic.developerLog(sb.toString());

    if (mEds != null) {
      for (final UnidirectionalEditDistance mEd : mEds) {
        mEd.logStats();
      }
    }

    if (LOG_AS_HISTO) {
      sb.setLength(0);
      sb.append("Alignment score histogram for all alignments attempted:").append(StringUtils.LS);
      for (int i = 0; i < mASHistogramUnderMaxScore.length - 2; i++) {
        if (mASHistogramUnderMaxScore[i] > 0 || mASHistogramOverMaxScore[i] > 0) {
          sb.append("AS:\t").append(i).append('\t').append(mASHistogramUnderMaxScore[i]).append('\t').append(mASHistogramOverMaxScore[i]).append(StringUtils.LS);
        }
      }
      if (mASHistogramOverMaxScore[mASHistogramOverMaxScore.length - 2] != 0) {
        sb.append("AS:\t255+\t0\t").append(mASHistogramOverMaxScore[mASHistogramOverMaxScore.length - 2]).append(StringUtils.LS);
      }
      if (mASHistogramUnderMaxScore[mASHistogramUnderMaxScore.length - 1] != 0) {
        sb.append("AS:\tMAXINT\t").append(mASHistogramUnderMaxScore[mASHistogramUnderMaxScore.length - 1]).append("\t0").append(StringUtils.LS);
      }
      Diagnostic.developerLog(sb.toString());
    }
  }
}
