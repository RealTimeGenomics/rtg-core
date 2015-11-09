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

import java.io.IOException;
import java.util.Arrays;
import java.util.Properties;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DNA;
import com.rtg.ngs.NgsParams;
import com.rtg.util.PropertiesUtils;
import com.rtg.util.QuickSort;
import com.rtg.util.QuickSortIntIntProxy;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Edit distance which aligns mismatches and at most one indel.
 */
@TestClass({"com.rtg.alignment.SingleIndelEditDistanceTest", "com.rtg.alignment.SingleIndelEditDistancePropsFileTest"})
public class SingleIndelEditDistance extends IntegralAbstract implements UnidirectionalEditDistance {

//  private static final String INDEL_TABLE = GlobalFlags.getStringValue(GlobalFlags.EDIT_DIST_INDEL_TABLE_FLAG);
  //private static final String INDEL_TABLE = System.getProperty(GlobalFlags.EDIT_DIST_INDEL_TABLE_FLAG);
  //private static final String INDEL_TABLE = "/home/len/indeltable.properties";

  protected static final int UNKNOWN = DNA.N.ordinal();

  protected final int mSubstitutionPenalty;
  protected final int mUnknownsPenalty;

  //The following variables are used temporarily during a call to calculateEditDistance
  //WARNING: this makes the code non-reentrant.
  private boolean mBestValid = false;
  //The following have not been initialized or set unless mBestValid == true
  protected byte[] mRead;
  protected byte[] mTemplate;
  protected int mRLen;
  protected int mMaxScore;
  protected int mMaxShift;
  protected final int[] mDiagonal;
  protected final int[] mDiagonalCum;
  protected int mFirstMiss;
  protected int mLastMiss;
  protected int mDiagScore;
  protected int mBestOffset;
  //The following are only set if bestOffset != 0
  protected int mBestPosn;
  protected int mBestScore;
  protected boolean mBestForward; //true iff the orientation is forward (the off diagonal offsets are searched from end backward)

  //used to return the actions array
  //WARNING: this makes the code non-reentrant.
  private int[] mWorkspace;

  int mMaxRLenMaxShift;

  private final int[] mIndelOffsets; // Table of insertion / deletion lengths to consider. Negative indicates deletion from reference
  private final int[] mIndelPenalties; // Alignment score penalty associated with the corresponding entry in mIndelOffsets.

    /**
     * Constructor
     * @param ngsParams containing configuration settings for alignment penalties
     * @param maxReadLength maximum read length in data set
     */
  public SingleIndelEditDistance(NgsParams ngsParams, int maxReadLength) {

    mIndelOffsets = new int[2 * maxReadLength];
    mIndelPenalties = new int[2 * maxReadLength];

    if (ngsParams.singleIndelPenalties() == null) {
      loadAffinePenalties(mIndelOffsets, mIndelPenalties, ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty());
      mSubstitutionPenalty = ngsParams.substitutionPenalty();
      mUnknownsPenalty = ngsParams.unknownsPenalty();
    } else {
      final Properties pr = loadIndelPenalties(mIndelOffsets, mIndelPenalties, ngsParams.singleIndelPenalties());
      try {
        mSubstitutionPenalty = Integer.parseInt(pr.getProperty("error_snp_penalty"));
        Diagnostic.developerLog("Loaded substitution penalty: " + mSubstitutionPenalty);
      } catch (RuntimeException e) {
        throw new NoTalkbackSlimException("Could not parse substitution penalty in penalties file " + ngsParams.singleIndelPenalties() + ": " + e.getMessage());
      }
      try {
        mUnknownsPenalty = Integer.parseInt(pr.getProperty("error_unknowns_penalty"));
        Diagnostic.developerLog("Loaded unknowns penalty: " + mUnknownsPenalty);
      } catch (RuntimeException e) {
        throw new NoTalkbackSlimException("Could not parse unknowns penalty in penalties file " + ngsParams.singleIndelPenalties() + ": " + e.getMessage());
      }
    }

    mDiagonal = new int[maxReadLength];
    mDiagonalCum = new int[maxReadLength];
  }

  /**
   * @param pos position on template
   * @return base value at given position, or 0 (N) if outside of template
   */
  byte template(int pos) {
    if (pos >= 0 && pos < mTemplate.length) {
      return mTemplate[pos];
    }
    return 0; // value of unknown.
  }

  static Properties loadIndelPenalties(int[] offsets, int[] penalties, String indelFile) {
    assert (offsets.length & 1) == 0;
    assert offsets.length == penalties.length;
    try {
      final Properties pr = PropertiesUtils.getPriorsResource(indelFile, PropertiesUtils.PropertyType.ALIGNMENT_PROPERTY_TYPE);
      // The property gives a comma separated list of integer penalties according to each length of insert, first element corresponding to insertion of length 1
      final int[] insdist = ArrayUtils.parseIntArray(pr.getProperty("error_ins_penalty"));
      if (insdist.length == 0) {
        throw new NumberFormatException("Invalid error_ins_penalty");
      }
      final int[] deldist = ArrayUtils.parseIntArray(pr.getProperty("error_del_penalty"));
      if (deldist.length == 0) {
        throw new NumberFormatException("Invalid error_del_penalty");
      }
      // These factors give the rate at which the penalty increases for each additional base past the last specifed offset penalty
      final double insslope = Double.valueOf(pr.getProperty("error_ins_penalty_extension_slope", "0.5"));
      final double delslope = Double.valueOf(pr.getProperty("error_del_penalty_extension_slope", "0.5"));
      int lastInsPenaltyOffset = 0;
      int lastDelPenaltyOffset = 0;
      for (int i = 0, offset = 0; i < offsets.length - 1; offset++) {
        offsets[i] = offset + 1;
        int penalty;
        if (offset < insdist.length) {
          penalty = insdist[offset];
          lastInsPenaltyOffset = offset;
        } else {
          penalty = insdist[lastInsPenaltyOffset] + (int) (insslope * (offset - lastInsPenaltyOffset)); // Extrapolate from last given penalty
        }
        penalties[i] = penalty;
        i++;

        offsets[i] = -offset - 1;
        if (offset < deldist.length) {
          penalty = deldist[offset];
          lastDelPenaltyOffset = offset;
        } else {
          penalty = deldist[lastDelPenaltyOffset] + (int) (delslope * (offset - lastDelPenaltyOffset)); // Extrapolate from last given penalty
        }
        penalties[i] = penalty;
        i++;
      }
      // Sort by increasing penalty.
      QuickSort.sort(new QuickSortIntIntProxy(penalties, offsets, true) {
        @Override
        public int compare(long index1, long index2) {
          final int res = super.compare(index1, index2);
          if (res != 0) {
            return res;
          }
          // Break ties by preferring smallest absolute offset, and if those are equal, prefer negative offset (deletion).
          final int o1 = Math.abs(mPairs[(int) index1]);
          final int o2 = Math.abs(mPairs[(int) index2]);
          return o1 != o2 ? Integer.compare(o1, o2) : Integer.compare(mPairs[(int) index1], mPairs[(int) index2]);
        }
      });
      Diagnostic.developerLog("Loaded indel table offsets: " + Arrays.toString(offsets));
      Diagnostic.developerLog("Loaded indel table penalties: " + Arrays.toString(penalties));
      return pr;
    } catch (NumberFormatException e) {
      throw new NoTalkbackSlimException("Could not parse entries in indel penalties file " + indelFile + ": " + e.getMessage());
    } catch (IOException e) {
      throw new NoTalkbackSlimException("Could not load indel penalties file " + indelFile + ": " + e.getMessage());
    }
  }

  static void loadAffinePenalties(int[] offsets, int[] penalties, int gapOpenPenalty, int gapExtendPenalty) {
    for (int i = 0; i < penalties.length; i++) {
      offsets[i] = i / 2 + 1;
      penalties[i] = gapOpenPenalty + offsets[i] * gapExtendPenalty;
      if ((i & 1) == 0) {
        offsets[i] = -offsets[i];
      }
    }
  }

  @Override
  public void logStats() {
    // do nothing
  }

  void resizeWorkspace(int size) {
    final int workspacesize = ActionsHelper.ACTIONS_START_INDEX + 1 + (int) (size / (double) ActionsHelper.ACTIONS_PER_INT + 0.5);
    mWorkspace = new int[workspacesize];
    mMaxRLenMaxShift = size;
  }

  /**
   * Tries to perform a quick alignment with substitutions only and as many
   * as possible, permitted by the penalties. This is not
   * useful for Complete Genomics data (and all Complete Genomics testing has
   * been removed from it).
   * <br><br>
   * The following diagram illustrates the relationship of the various parts.
   * More detailed diagrams are given for individual methods below.
   * <br>
   * <img src="doc-files/singleIndelGeneral.jpg" alt="explanatory image">
   * @param read the encoded read bases
   * @param rLen length of read (assumed that given read is at least this long)
   * @param template the encoded template
   * @param zeroBasedStart start position in the template
   * @param maxScore helps with early termination.
   * @param maxShift maximum amount that start position can be shifted by (also, the length of the longest indel that can be found).
   * @param cgLeft ignored.
   * @return the actions array, or null if not possible
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rLen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
//    System.err.println("rLen=" + rLen + " template.length=" + template.length + " zeroBasedStart=" + zeroBasedStart + " maxScore=" + maxScore + " maxShift=" + maxShift);
    assert rLen <= mDiagonal.length;
    mRead = read;
    mTemplate = template;

    if (rLen + maxShift > mMaxRLenMaxShift) {
      resizeWorkspace(rLen + maxShift);
    }

    mRLen = rLen;
    mMaxScore = maxScore;
    mMaxShift = maxShift;

    initDiagonal(rLen, zeroBasedStart);

    boolean insBailout = false;
    boolean delBailout = false;
    for (int i = 0; i < mIndelOffsets.length && !(insBailout && delBailout); i++) {
      final int offset = mIndelOffsets[i];

      if (offset > 0 && offset > maxShift) {
        insBailout = true;
        continue;
      } else if (offset < 0 && -offset > maxShift) {
        delBailout = true;
        continue;
      }

      final int offScore = mIndelPenalties[i];
//      System.err.println("mBestScore=" + mBestScore + " mBestOffset=" + mBestOffset + " mBestPosn=" + mBestPosn + " mBestOrientation=" + mBestForward);
      if (offScore >= mBestScore) {
        break;
      }

      if (offset > 0) {
        final int firsti = mFirstMiss + offset;
        final int tEndFN = zeroBasedStart - offset + rLen;
        if (firsti <= rLen) {
          offDiagonalForwardNegative(firsti, rLen, tEndFN, -offset, offScore, mDiagScore - mDiagonalCum[rLen - offset]);
          if (offScore >= mBestScore) {
            break;
          }
        }

        final int rEnd = Math.min(mLastMiss + 1, rLen - offset);
        //System.err.println("i=" + i + " rLength=" + read.length + " rLen=" + rLen + " tEndRP=" + tEndRP + " tlength=" + template.length);
        if (rEnd > 0) {
          offDiagonalReversePositive(0, rEnd, zeroBasedStart + offset, offset, offScore);
          if (offScore >= mBestScore) {
            break;
          }
        }
      } else {

        final int tEndFP = zeroBasedStart - offset + rLen;
          offDiagonalForwardPositive(mFirstMiss, rLen, tEndFP, -offset, offScore);
          if (offScore >= mBestScore) {
            break;
          }
        final int tStart = zeroBasedStart + offset;
        offDiagonalReverseNegative(0, mLastMiss + 1, tStart, offset, offScore);
      }
    }
//    System.err.println(toString());
    //System.err.println("mBestScore=" + mBestScore + " mBestOffset=" + mBestOffset + " mBestPosn=" + mBestPosn + " mBestOrientation=" + mBestForward);
    final int[] actions = actions(zeroBasedStart, rLen);
    if (actions != null) {
      actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = mBestScore;
    }
    assert actions == null || ActionsHelper.alignmentScore(actions) == mBestScore;
    return actions;
  }

  protected void initDiagonal(int rLen, int zeroBasedStart) {
    mDiagScore = 0;
    int fMis = -1; //position of first mismatch
    int lMis = -1; //position of last mismatch
    for (int i = 0; i < rLen; i++) {
      mDiagonalCum[i] = mDiagScore;
      final byte tb = template(zeroBasedStart + i); //mTemplate[zeroBasedStart + i];
      final byte rb = mRead[i];
      final int m = m(rb, tb, mSubstitutionPenalty, mUnknownsPenalty);
      mDiagonal[i] = m;
      mDiagScore += m;
      if (m > 0) {
        if (fMis == -1) {
          fMis = i;
        }
        lMis = i;
      }
    }
    mFirstMiss = fMis < 0 ? 0 : fMis;
    mLastMiss = lMis < 0 ? rLen - 1 : lMis;

    mBestOffset = 0;
    mBestScore = Math.min(mDiagScore, mMaxScore + (mMaxScore == Integer.MAX_VALUE ? 0 : 1)); //have to BEAT best score, but our maxscore is an inclusive value, hence +1
    mBestValid = true;
  }

  /**
   * @param zeroBasedStart start on template.
   * @param rLen read length
   * @return the actions array using the best values.
   */
  protected int[] actions(int zeroBasedStart, int rLen) {
    mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = 0;
    mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
    if (mBestOffset == 0) { //potentially on the main diagonal
      if (mBestScore <= mMaxScore) {
        mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart + rLen; //add rlen because the actionshelper prepend "helpfully" alters the workspace start position as it goes
        //use the diagonal
        action(0, rLen, zeroBasedStart + rLen);
        return mWorkspace;
      } else {
        //nothing bettered maxscore
        return null;
      }
    }

    assert mBestOffset != 0;
    if (mBestForward) { //Forward
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart + rLen + mBestOffset; //add rlen because the actionshelper prepend "helpfully" alters the workspace start position as it goes
      action(mBestPosn, rLen, zeroBasedStart + rLen + mBestOffset);
      indelAction();
      final int rEnd = mBestPosn + (mBestOffset > 0 ? 0 : mBestOffset);
      action(0, rEnd, zeroBasedStart + rEnd);
    } else { //reverse
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart + rLen; //add rlen because the actionshelper prepend "helpfully" alters the workspace start position as it goes
      final int rStart = mBestPosn + (mBestOffset > 0 ? mBestOffset : 0) + 1;
      action(rStart, rLen, zeroBasedStart + rLen);
      indelAction();
      action(0, mBestPosn + 1, zeroBasedStart + mBestOffset + mBestPosn + 1);
    }
    return mWorkspace;
  }

  private void indelAction() {
    assert mBestOffset != 0;
    if (mBestOffset > 0) {
      indelAction(mBestOffset, mBestForward ? ActionsHelper.DELETION_FROM_REFERENCE : ActionsHelper.INSERTION_INTO_REFERENCE);
    } else {
      indelAction(-mBestOffset, mBestForward ? ActionsHelper.INSERTION_INTO_REFERENCE : ActionsHelper.DELETION_FROM_REFERENCE);
    }
  }

  private void indelAction(final int cnt, final int act) {
    //System.err.println("indelAction cnt=" + cnt + " act=" + act + " score=" + score);
    ActionsHelper.prepend(mWorkspace, cnt, act, 0);
  }

  /**
   * @param rStart start position on read (zero based, inclusive)
   * @param rEnd end position on read (zero based, exclusive)
   * @param tEnd  end position on template (zero based, exclusive)
   */
  private void action(final int rStart, final int rEnd, final int tEnd) {
    //System.err.println("action rStart=" + rStart + " rEnd=" + rEnd + " tEnd=" + tEnd);
    if (rStart == rEnd) {
      return;
    }
    int cnt = 0;
    int lastScore = 0;
    int lastAction = -1;
    for (int r = rEnd - 1, t = tEnd - 1; r >= rStart; r--, t--) {

      final byte tNt = template(t);
      final byte rNt = mRead[r];

      final int thisAction;
      final int thisScore;

      //Do some ridiculous crap to work out if the last action was the same, mismatch, or unknown.
      if (tNt * rNt == UNKNOWN) {
        thisAction = ActionsHelper.MISMATCH;
        thisScore = mUnknownsPenalty;
      } else {
        if (tNt == rNt) {
          thisAction = ActionsHelper.SAME;
          thisScore = 0;
        } else {
          thisAction = ActionsHelper.MISMATCH;
          thisScore = mSubstitutionPenalty;
        }
      }

      if (thisScore == lastScore && thisAction == lastAction) {
        cnt++;
      } else {
        prepend(cnt, lastScore, lastAction); //first time through, cnt == 0 which will be skipped over
        lastScore = thisScore;
        lastAction = thisAction;
        cnt = 1;
      }
    }
    prepend(cnt, lastScore, lastAction);
  }

  private void prepend(final int cnt, final int actionPenalty, final int action) {
    if (cnt > 0) {
      //System.err.println("prepend cnt=" + cnt + " lastM=" + lastM + " act=" + act);
      ActionsHelper.prepend(mWorkspace, cnt, action, cnt * actionPenalty);
    }
  }

  /**
   * <img src="doc-files/singleIndelForwardPositive.jpg" alt="explanatory image">
   * @param rStart start position for part of read to be explored (zero based).
   * @param rEnd end position for part of read to be explored (zero based, exclusive).
   * @param tEnd  end position for template (zero based, exclusive).
   * @param offSet distance to the diagonal being explored (&ne; 0).
   * @param offsetPenalty the alignment score penalty associated with this offset
   */
  protected void offDiagonalForwardPositive(int rStart, int rEnd, int tEnd, int offSet, int offsetPenalty) {
    assert rStart < rEnd;
    assert offSet > 0;
    //System.err.println("offDiagonalForwardPositive: rStart=" + rStart + " rEnd=" + rEnd + " tEnd=" + tEnd + " offset=" + offSet);
    assert offSet != 0;
    int score = mDiagScore + offsetPenalty;
    int possScore = offsetPenalty;
    for (int r = rEnd - 1, t = tEnd - 1; r >= rStart; r--, t--) {
      final int m = match(r, t);
      score += m - mDiagonal[r];
      possScore += m;
      if (possScore >= mBestScore) {
        break;
      }
      if (score < mBestScore) {
        mBestScore = score;
        mBestOffset = offSet;
        mBestPosn = r;
        mBestForward = true;
      }
    }
  }

  /**
   * <img src="doc-files/singleIndelForwardNegative.jpg" alt="explanatory image" >
   * @param rStart start position for read (zero based).
   * @param rEnd end position for read (zero based, exclusive).
   * @param tEnd  end position for template (zero based, exclusive).
   * @param offSet distance to the diagonal being explored (&ne; 0).
   * @param offsetPenalty the alignment score penalty associated with this offset
   * @param scoreAdjust adjustment to the score when there are deletions on the read.
   */
  protected void offDiagonalForwardNegative(int rStart, int rEnd, int tEnd, int offSet, int offsetPenalty, int scoreAdjust) {
//    System.err.println("offDiagonalForwardNegative: rStart=" + rStart + " rEnd=" + rEnd + " tEnd=" + tEnd + " offset=" + offSet + " offsetPenalty=" + offsetPenalty + " scoreAdjust=" + scoreAdjust);
    assert rStart <= rEnd;
    assert offSet < 0;
    assert scoreAdjust >= 0;
    int score = mDiagScore + offsetPenalty - scoreAdjust;
    int possScore = offsetPenalty;

    if (score < mBestScore) {
      mBestScore = score;
      mBestOffset = offSet;
      mBestPosn = rEnd;
      mBestForward = true;
    }

//    System.err.println("oDFN score=" + score + " possScore=" + possScore);
    for (int r = rEnd - 1, t = tEnd - 1; r >= rStart; r--, t--) {
      final int m = match(r, t);
      final int diag = mDiagonal[r + offSet];
      score += m - diag;
      possScore += m;
//      System.err.println("oDFN " + " r=" + r + " m=" + m + " diag=" + diag + " score=" + score + " possScore=" + possScore + " mBestScore=" + mBestScore);
      if (possScore >= mBestScore) {
        break;
      }
      if (score < mBestScore) {
        mBestScore = score;
        mBestOffset = offSet;
        mBestPosn = r;
        mBestForward = true;
      }
    }
  }

  /**
   * <img src="doc-files/singleIndelReversePositive.jpg" alt="explanatory image" >
   * @param rStart start position for read (zero based).
   * @param rEnd end position to be explored on read (zero based, exclusive).
   * @param tStart  start position for template (zero based).
   * @param offSet distance to the diagonal being explored (&gt; 0).
   * @param offsetPenalty the alignment score penalty associated with this offset
   */
  protected void offDiagonalReversePositive(int rStart, int rEnd, int tStart, int offSet, int offsetPenalty) {
    assert rStart < rEnd;
    assert offSet > 0;
    int score = mDiagScore + offsetPenalty - mDiagonalCum[offSet];
    int possScore = offsetPenalty;

//    System.err.println("offDiagonalReversePositive: rStart=" + rStart + " rEnd=" + rEnd + " tStart=" + tStart + " offset=" + offSet + " score=" + score + " possScore=" + possScore + " mDiagScore=" + mDiagScore);

    if (score < mBestScore) {
      mBestScore = score;
      mBestOffset = offSet;
      mBestPosn = -1;
      mBestForward = false;
    }

    for (int r = rStart, t = tStart; r < rEnd; r++, t++) {
      final int m = match(r, t);
      final int diag = mDiagonal[r + offSet];
      score += m - diag;
      possScore += m;
//      System.err.println("oDRP " + " r=" + r + " t=" + t + " m=" + m + " diag=" + diag + " score=" + score + " possScore=" + possScore + " mBestScore" + mBestScore);
      if (possScore >= mBestScore) {
        break;
      }
      if (score < mBestScore) {
        mBestScore = score;
        mBestOffset = offSet;
        mBestPosn = r;
        mBestForward = false;
      }
    }
  }

  /**
   * <img src="doc-files/singleIndelReverseNegative.jpg" alt="explanatory image" >
   * @param rStart start position for read (zero based).
   * @param rEnd end position for read (zero based, exclusive).
   * @param tStart  start position for template (zero based).
   * @param offSet distance to the diagonal being explore (&ne; 0).
   * @param offsetPenalty the alignment score penalty associated with this offset
   */
  protected void offDiagonalReverseNegative(int rStart, int rEnd, int tStart, int offSet, int offsetPenalty) {
    //System.err.println("offDiagonalReverseNegative: rStart=" + rStart + " rEnd=" + rEnd + " tStart=" + tStart + " offset=" + offSet);
    assert rStart <= rEnd;
    assert offSet < 0;
    int score = mDiagScore + offsetPenalty;
    int possScore = offsetPenalty;
    for (int r = rStart, t = tStart; r < rEnd; r++, t++) {
      final int m = match(r, t);
      final int diag = mDiagonal[r];
      score += m - diag;
      possScore += m;
      //System.err.println("oDRN " + " r=" + r + " t=" + t + " m=" + m + " diag=" + diag + " score=" + score + " possScore=" + possScore);
      if (possScore >= mBestScore) {
        break;
      }
      if (score < mBestScore) {
        mBestScore = score;
        mBestOffset = offSet;
        mBestPosn = r;
        mBestForward = false;
      }
    }
  }

  /**
   * @param r position in read (zero based).
   * @param t position in template (zero based).
   * @return score for a match or mismatch between this pair.
   */
  private int match(int r, int t) {
    //System.err.println("match: r=" + r + " t=" + t);
    final byte tb = template(t); //mTemplate[t];
    final byte rb = mRead[r];
    final int m;
    m = m(tb, rb, mSubstitutionPenalty, mUnknownsPenalty);
    //System.err.println("match: r=" + r + " t=" + t + " tb=" + tb + " rb=" + rb + " m=" + m);
    return m;
  }

  static int m(final byte tb, final byte rb, final int mismatchPenalty, final int unknownsPenalty) {
    assert UNKNOWN == 0;
    if (rb * tb == UNKNOWN) {
      return unknownsPenalty;
    }
    return rb == tb ? 0 : mismatchPenalty;
  }

  /**
   * @return copy of internal indel penalties. For testing purposes.
   */
  protected int[] getIndelPenalties() {
    return Arrays.copyOf(mIndelPenalties, mIndelPenalties.length);
  }
  /**
   * @return copy of internal indel penalties. For testing purposes.
   */
  protected int[] getIndelOffsets() {
    return Arrays.copyOf(mIndelOffsets, mIndelOffsets.length);
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

  @Override
  //This is intended mainly for debugging not for external display
  public void toString(StringBuilder sb) {
    sb.append("SingleIndelEditDistance");
    sb.append(" unknownsPenalty=").append(mUnknownsPenalty);
    sb.append(" substitution=").append(mSubstitutionPenalty);
    sb.append(LS);

    if (!mBestValid) {
      return;
    }
    sb.append(" rLen=").append(mRLen);
    sb.append(" maxScore=").append(mMaxScore);
    sb.append(" mMaxShift=").append(mMaxShift);
    sb.append(LS);

    for (int i = 0; i < mRLen; i++) {
      sb.append(" ").append(mDiagonal[i]);
    }
    sb.append(LS);

    for (int i = 0; i < mRLen; i++) {
      sb.append(" ").append(mDiagonalCum[i]);
    }
    sb.append(LS);

    sb.append(" firstMiss=").append(mFirstMiss);
    sb.append(" lastMiss=").append(mLastMiss);
    sb.append(" diagScore=").append(mDiagScore);
    sb.append(" bestOffset=").append(mBestOffset);
    sb.append(LS);

    if (mBestOffset == 0) {
      return;
    }
    sb.append(" bestPosn=").append(mBestPosn);
    sb.append(" bestScore=").append(mBestScore);
    sb.append(" bestForward=").append(mBestForward);
    sb.append(LS);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    int tot = 0;
    for (int i = 0; i < mRLen; i++) {
      final int diag = mDiagonal[i];
      Exam.assertTrue(diag == 0 || diag == mSubstitutionPenalty);
      Exam.assertEquals(tot, mDiagonalCum[i]);
      tot += diag;
    }
    Exam.assertEquals(tot, mDiagScore);
    return true;
  }


  @Override
  public boolean integrity() {
    Exam.assertEquals(0, UNKNOWN);
    Exam.assertTrue(mSubstitutionPenalty > 0);
    Exam.assertNotNull(mIndelOffsets);
    Exam.assertNotNull(mIndelPenalties);
    Exam.assertEquals(mIndelOffsets.length, mIndelPenalties.length);
    int lastPenalty = 0;
    int bestPenalty = 0;
    for (int i = 0; i < mIndelPenalties.length; i++) {
      final int current = mIndelPenalties[i];
      Exam.assertTrue(current >= lastPenalty);
      lastPenalty = current;
      if (mIndelOffsets[i] == mBestOffset) {
        bestPenalty = mIndelPenalties[i];
      }
    }
    Exam.assertNotNull(mDiagonal);
    Exam.assertNotNull(mDiagonalCum);
    if (!mBestValid) {
      return true;
    }
    Exam.assertNotNull(mRead);
    Exam.assertNotNull(mTemplate);
    Exam.assertTrue(0 <= mRLen && mRLen <= mDiagonal.length);
    Exam.assertTrue(mMaxScore >= 0);
    Exam.assertTrue(mMaxShift >= 0);
    Exam.assertTrue(0 <= mFirstMiss && mFirstMiss < mRLen);
    if (mFirstMiss > 0) {
      Exam.assertEquals(mSubstitutionPenalty, mDiagonal[mFirstMiss]);
    }
    Exam.assertTrue(0 <= mLastMiss && mLastMiss < mRLen);
    if (mLastMiss < mRLen - 1) {
      Exam.assertEquals(mSubstitutionPenalty, mDiagonal[mLastMiss]);
    }
    Exam.assertTrue(0 <= mDiagScore && mDiagonalCum[mRLen - 1] <= mDiagScore);
    if (mBestOffset == 0) {
      return true;
    }
    Exam.assertTrue(0 <= mBestPosn && mBestPosn <= mRLen);
    Exam.assertTrue(bestPenalty <= mBestScore);
    Exam.assertTrue(mBestScore <= mDiagScore && mBestScore <= mMaxScore);

    return true;
  }

}
