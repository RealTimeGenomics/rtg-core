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

import java.util.Arrays;

import com.rtg.ngs.NgsParams;
import com.rtg.util.integrity.Exam;

/**
 * Edit distance which aligns mismatches and at most one indel.
 * Uses seeding to achieve O(N).
 */
public class SingleIndelSeededEditDistance extends SingleIndelEditDistance {

  private static final int SEED_COUNT = 0;
  private static final int SEED_LAST = 1;
  private static final int SEED_LENGTH = SEED_LAST + 1;

  private final boolean mFindDiagonal;
  private final int mDeltaThreshold;
  private final int mDiffThreshold;

  private final Seed mSeed;

  /** Counts and last seen pointers for seeds. */
  private final short[] mSeedInfo;

  private int[] mDelta = new int[21]; //long enough for common case of 100 long reads and  max shift of 10
  private int mDeltaLength;

  private final int[] mCountAmbiguous = new int[5];
  private int mAmbiguousRescued = 0;

  protected final int mIndelOpenPenalty;

  protected final int mIndelExtendPenalty;


  /**
   * Default parameters.
   * @param ngsParams parameters as supplied in command line etc.
   * @param maxReadLength maximum read length in data set
   */
  public SingleIndelSeededEditDistance(final NgsParams ngsParams, int maxReadLength) {
    this(ngsParams, 3, 2, 2, false, maxReadLength);
  }

  /**
   * @param ngsParams parameters as supplied in command line etc.
   * @param seedSize size of seeds. Size is heuristic and probably should be in the range 3 to 5. Possibly should increase with (log of?) increasing read length.
   * @param deltaThreshold the minimum number of seed hits on an off main diagonal before it will be considered. Increasing this will reduce the number of alignments found
   *        (which may be bad) and also reduce time taken.
   * @param diffThreshold the minimum difference in mismatch count between upper half and lower half of read before a direction will be chosen and off diagonal processing done.
   * @param findDiagonal if set, do not assume there is a good hit along the diagonal, do an initial search for the best diagonal. Then pass alignment position to single indel edit distance
   * @param maxReadLength maximum read length in data set
   */
  public SingleIndelSeededEditDistance(final NgsParams ngsParams, final int seedSize, final int deltaThreshold, final int diffThreshold, boolean findDiagonal, int maxReadLength) {
    super(ngsParams, maxReadLength);
    mSeed = new Seed(seedSize);
    mSeedInfo = new short[SEED_LENGTH * (1 << mSeed.size() * 2)];
    mDelta = new int[0];
    mDeltaThreshold = deltaThreshold;
    mDiffThreshold = diffThreshold * mSubstitutionPenalty;
    mIndelOpenPenalty = ngsParams.gapOpenPenalty();
    mIndelExtendPenalty = ngsParams.gapExtendPenalty();
    mFindDiagonal = findDiagonal;
  }

  private SingleIndelEditDistance mSied = null;
  private GotohEditDistance mGotoh = null;

  /**
   * Forces searching for best diagonal, and then passing off that to the single indel edit distance.
   * @param ngsParams parameters as supplied in command line etc.
   * @param useGotoh true to use Gotoh edit distance to attempt to retain alignments that the single indel edit distance failed to
   * @param seedSize size of seeds. Size is heuristic and probably should be in the range 3 to 5. Possibly should increase with (log of?) increasing read length.
   * @param deltaThreshold the minimum number of seed hits on an off main diagonal before it will be considered. Increasing this will reduce the number of alignments found
*        (which may be bad) and also reduce time taken.
   * @param diffThreshold the minimum difference in mismatch count between upper half and lower half of read before a direction will be chosen and off diagonal processing done.
   * @param maxReadLength maximum read length in data set
   */
  public SingleIndelSeededEditDistance(NgsParams ngsParams, boolean useGotoh, int seedSize, int deltaThreshold, int diffThreshold, int maxReadLength) {
    this(ngsParams, seedSize, deltaThreshold, diffThreshold, true, maxReadLength);
    mSied = new SingleIndelEditDistance(ngsParams, maxReadLength);
    mGotoh = useGotoh ? new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false) : null;
  }

  /**
   * The following diagram illustrates the relationship of the various parts.
   * More detailed diagrams are given for individual methods below.
   * <br>
   * <img src="doc-files/singleIndelGeneral.jpg" alt="explanatory image">
   * @param read the encoded read bases
   * @param rLen length of read (assumed that given read is at least this long)
   * @param template the encoded template
   * @param initialZeroBasedStart start position in the template
   * @param maxScore helps with early termination.
   * @param maxShift ignored. No shifting permitted.
   * @param cgLeft ignored.
   * @return the actions array, or null if not possible
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rLen, byte[] template, int initialZeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    //System.err.println("rLen=" + rLen + " template.length=" + template.length + " zeroBasedStart=" + zeroBasedStart + " maxScore=" + maxScore + " maxShift=" + maxShift);
    assert rLen <= mDiagonal.length;
    // Quick idiocy checks first to avoid string construction cost
//    if (initialZeroBasedStart < 0 || initialZeroBasedStart + rLen > template.length) {    //if read mapped off template at either end, just delegate to next in chain
//      return null;
//    }
    mRead = read;
    mTemplate = template;

    if (rLen + maxShift > mMaxRLenMaxShift) {
      resizeWorkspace(rLen + maxShift);
    }

    mRLen = rLen;
    mMaxScore = maxScore;
    mMaxShift = maxShift;
    mDeltaLength = 2 * mMaxShift + 1;
    if (mDelta.length < (mDeltaLength + 2)) {
      mDelta = new int[mDeltaLength + 2];
    }
    int zeroBasedStart = initialZeroBasedStart;

    initDiagonal(rLen, zeroBasedStart);

    if (mDiagScore <= mIndelOpenPenalty + mIndelExtendPenalty) {
      //no point in looking for off diagonals - go for the diagonal
      return actions(zeroBasedStart, rLen);
    }

    populateDeltas(read, rLen, template, zeroBasedStart);

//    System.err.println(" d " + deltaToString());

    if (mFindDiagonal) {
      //System.err.println(deltaToString());
      //Look for a unique maximal offset
      int deltaOff = 0;
      int deltaBest = 0;
      int deltaCnt = 0;
      int deltaOtherOff = 0;
      for (int i = 0; i <= mDeltaLength + 1; ++i) {
        final int delta = mDelta[i];
        if (delta > deltaBest) {
          deltaCnt = 1;
          deltaBest = delta;
          deltaOff = i - maxShift - 1;
//          System.err.println("new best delta " + delta + " @ " + deltaOff);
        } else if (delta == deltaBest) {
          ++deltaCnt;
          deltaOtherOff = i - maxShift - 1;
//          System.err.println("found equal best delta " + delta + " @ " + (i - maxShift - 1));
//        } else if (delta > 10) {
//          System.err.println("delta " + delta + " @ " + (i - maxShift - 1));
        }
      }

      if (deltaBest != 0) {
        if (mSied != null) {
          if (deltaCnt == 0) {
            return null;
          } else if (deltaCnt == 1) {
            return hybridAlignment(read, rLen, template, zeroBasedStart + deltaOff, maxScore, maxShift, cgLeft);
          } else if (deltaCnt == 2) {
            mCountAmbiguous[0]++;

            //using the gotoh aligner as a backup in this case doesn't seem to buy us a lot.

            int[] actions = mSied.calculateEditDistance(read, rLen, template, zeroBasedStart + deltaOff, maxScore, maxShift, cgLeft);
//            int[] actions = hybridAlignment(read, rLen, template, zeroBasedStart + deltaOff, maxScore, maxShift, cgLeft);
            final int[] firstActions = actions == null ? null : Arrays.copyOf(actions, actions.length);
            actions = mSied.calculateEditDistance(read, rLen, template, zeroBasedStart + deltaOtherOff, maxScore, maxShift, cgLeft);
//            actions = hybridAlignment(read, rLen, template, zeroBasedStart + deltaOtherOff, maxScore, maxShift, cgLeft);
            return tiebreakAmbiguous(firstActions, actions, maxScore);

          } else {
            if (deltaCnt > 5) {
              mCountAmbiguous[4]++;
            } else {
              mCountAmbiguous[deltaCnt - 2]++;
            }
          }
          return null;
        }


        zeroBasedStart += deltaOff;
//        System.err.println("zerbs " + zeroBasedStart);

        if (zeroBasedStart < 0 || zeroBasedStart + rLen > template.length) {    //if read mapped off template at either end, just delegate to next in chain
          return null;
        }
        initDiagonal(rLen, zeroBasedStart);
        if (mDiagScore <= mIndelOpenPenalty + mIndelExtendPenalty) {
          //no point in looking for off diagonals - go for the diagonal
          return actions(zeroBasedStart, rLen);
        }
        populateDeltas(read, rLen, template, zeroBasedStart);
      }
      return null;
    }


//    System.err.println(deltaToString());
    //Look for a unique maximal offset
    int deltaOff = 0;
    int deltaCnt = 0;
    int deltaBest = 0;
    for (int i = 0; i <= mDeltaLength + 1; ++i) {
      if (i == (maxShift + 1)) { //skip zero offset
        continue;
      }
      final int delta = mDelta[i];
      if (delta > deltaBest) {
        deltaBest = delta;
        deltaCnt = 1;
        deltaOff = i - maxShift - 1;
//        System.err.println("current best pos: " + deltaOff);
      } else if (delta == deltaBest) {
//        System.err.println("equal best pos: " + deltaOff);
        ++deltaCnt;
      }
    }
//    System.err.println("deltaOff=" + deltaOff + " deltaCnt=" + deltaCnt + " deltaBest=" + deltaBest);
    if (deltaOff != 0 && deltaBest > mDeltaThreshold) {
      if (deltaCnt == 1) {
        doCall(zeroBasedStart, rLen, deltaOff, template.length);
      } else if (deltaCnt > 1) {
        if (deltaCnt > 5) {
          mCountAmbiguous[4]++;
        } else {
          mCountAmbiguous[deltaCnt - 2]++;
        }
      }
    }

    //System.err.println("mBestScore=" + mBestScore + " mBestOffset=" + mBestOffset + " mBestPosn=" + mBestPosn + " mBestOrientation=" + mBestForward);
    final int[] actions = actions(zeroBasedStart, rLen);
    if (actions != null) {
      actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = mBestScore;
    }
    return actions;
  }


  private int[] tiebreakAmbiguous(int[] actions1, int[] actions2, int maxScore) {
    if (actions1 == null) {
      if (actions2 != null && ActionsHelper.alignmentScore(actions2) < maxScore) {
        ++mAmbiguousRescued;
      }
      return actions2; //if both are null, this will return null.
    } else if (actions2 == null) {
      if (ActionsHelper.alignmentScore(actions1) < maxScore) {
        ++mAmbiguousRescued;
      }
      return actions1;
    }
    if (ActionsHelper.alignmentScore(actions2) == ActionsHelper.alignmentScore(actions1)) {
      return null;
    }
    if (ActionsHelper.alignmentScore(actions2) < ActionsHelper.alignmentScore(actions1)) {
      if (ActionsHelper.alignmentScore(actions2) < maxScore) {
        ++mAmbiguousRescued;
      }
      return actions2;
    }
    if (ActionsHelper.alignmentScore(actions1) < maxScore) {
      ++mAmbiguousRescued;
    }
    return actions1;
  }


  private int[] hybridAlignment(byte[] read, int rLen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
//    System.err.println("Attempting hybrid alignment @ " + zeroBasedStart);
    int[] actions = mSied.calculateEditDistance(read, rLen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
    if (mGotoh != null && actions == null) {
      actions = mGotoh.calculateEditDistance(read, rLen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
//      System.err.println("Gotoh aligned read " + DnaUtils.bytesToSequenceIncCG(read) + " from " + zeroBasedStart + " as " + ActionsHelper.toString(actions) + " at " + ActionsHelper.zeroBasedTemplateStart(actions));
    }
    return actions;
  }

  private void doCall(int zeroBasedStart, int rLen, int deltaOff, int templateLength) {
    final int half = mDiagonalCum[mRLen / 2];
    final int diff = mDiagScore - 2 * half;
    final int dir; //the logic below is heuristic
    if (diff < -mDiffThreshold) {
      dir = -1;
    } else {
      if (diff > mDiffThreshold) {
        dir = +1;
      } else {
        dir = 0;
      }
    }
    if (dir != 0) { //TODO heuristic - should it scale with read length?
      final int offsetScore = mIndelOpenPenalty + mIndelExtendPenalty * Math.abs(deltaOff);
      if (dir == +1) {
        if (deltaOff > 0) {
          final int tEndFP = zeroBasedStart + deltaOff + rLen;
          //System.err.println("tEndFP=" + tEndFP + " tlength=" + template.length);
          if (tEndFP <= templateLength) {
            offDiagonalForwardPositive(mFirstMiss, rLen, tEndFP, deltaOff, offsetScore);
          }
        } else {
          assert deltaOff < 0;
          final int firstI = mFirstMiss - deltaOff;
          final int tEndFN = zeroBasedStart + deltaOff + rLen; //check if this should be + deltaOff or - deltaOff
          if (firstI < rLen && tEndFN <= templateLength) {
            offDiagonalForwardNegative(firstI, rLen, tEndFN, deltaOff, offsetScore, mDiagonalCum[rLen + deltaOff]);
          }

        }
      } else {
        final int tEndRP = rLen + zeroBasedStart;
        assert dir == -1;
        final int rEnd = Math.min(mLastMiss + 1, rLen - deltaOff);
        if (tEndRP <= templateLength && deltaOff > 0 && rEnd > 0) {
          offDiagonalReversePositive(rEnd, zeroBasedStart + deltaOff, deltaOff, offsetScore);
        } else {
          final int tStart = zeroBasedStart + deltaOff;
          final int tEndRN = zeroBasedStart + rLen + deltaOff;
          if (tStart >= 0 && tEndRN <= templateLength) {
            offDiagonalReverseNegative(mLastMiss + 1, tStart, deltaOff, offsetScore);
          }
        }
      }
    }
  }

  private void populateDeltas(byte[] read, int rLen, byte[] template, int zeroBasedStart) {
    Arrays.fill(mDelta, 0);
    Arrays.fill(mSeedInfo, (short) 0);
    //initialize seeds
    final SeedShifter rSeed = new SeedShifter(mSeed, read, rLen, 0);
    final SeedShifter sSeed = new SeedShifter(mSeed, read, rLen, -mDeltaLength - 1);
    final SeedShifter tSeed = new SeedShifter(mSeed, template, template.length, zeroBasedStart - mMaxShift - 1);

    for (int i = 0; i <= rLen + mDeltaLength; ++i) {
      final int r = rSeed.next();
      if (rSeed.isValid()) {
        final int ix = SEED_LENGTH * r;
        mSeedInfo[ix + SEED_COUNT]++;
        mSeedInfo[ix + SEED_LAST] = (short) rSeed.position();
      }

      final int s = sSeed.next();
      if (sSeed.isValid()) {
        final int ix = SEED_LENGTH * s;
        final int cnt = mSeedInfo[ix + SEED_COUNT]--;
        assert cnt >= 0;
      }

      final int t = tSeed.next();
      if (tSeed.isValid()) {
        final int ix = SEED_LENGTH * t;
        final int cnt = mSeedInfo[ix + SEED_COUNT];
        if (cnt == 1) { //single unique hit in the read
          final int rp = mSeedInfo[ix + SEED_LAST];
          final int delta = tSeed.position() - zeroBasedStart - rp + mMaxShift + 1; //??offset
          mDelta[delta]++;
        }
      }
    }
  }

  @Override
  public boolean globalIntegrity() {
    super.globalIntegrity();
    for (int i = 0; i < mSeedInfo.length / SEED_LENGTH; ++i) {
      Exam.assertEquals(0, mSeedInfo[SEED_LENGTH * i + SEED_COUNT]);
    }
    for (int aMDelta : mDelta) {
      Exam.assertTrue(aMDelta >= 0);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertNotNull(mSeed);
    Exam.assertTrue(mDeltaThreshold >= 0);
    Exam.assertTrue(mDiffThreshold >= 0);
    return true;
  }


  @Override
  public void logStats() {
    super.logStats();
    System.out.println("SingleIndelSeededAligner skipped due to ambiguity: ");
    System.out.println("Count 2: " + mCountAmbiguous[0] + ", " + mAmbiguousRescued + " rescued");
    System.out.println("Count 3:  " + mCountAmbiguous[1]);
    System.out.println("Count 4:  " + mCountAmbiguous[2]);
    System.out.println("Count 5:  " + mCountAmbiguous[3]);
    System.out.println("Count 6+: " + mCountAmbiguous[4]);
  }
}
