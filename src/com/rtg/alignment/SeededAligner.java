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

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.ngs.NgsParams;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Fast seeded aligner, using incremental hashing, hash tables, bounding and
 * soft acceptance.
 *
 */
final class SeededAligner implements UnidirectionalEditDistance {

//  private static final Integer CHECK_ALIGNMENTS_EVERY_X = GlobalFlags.getIntegerValue(GlobalFlags.SEEDED_ALIGNER_CHECK_ALIGNMENTS_EVERY_X_FLAG);

  //warning: setting the below to true fixes some bad alignments, but drastically reduces the number of alignments seeded aligner "eats"
  private final boolean mNullOnLowComplexity;

  private static final int SEEDTYPE = 1;
  private static final int GAPTYPE = 0;

  static final int NTHASHSIZE = 15; // if this is over 15 then you'll be stuffed!
  private static final long MAXHASHVALUE = (int) Math.pow(4, NTHASHSIZE) - 1;

  private int mHashtableSize;
  private long[] mIndex;
  private int[] mPos;

  private final UnidirectionalEditDistance mED;
  private final UnidirectionalEditDistance mFixedStart;
  private final UnidirectionalEditDistance mFixedEnd;
//  private final UnidirectionalEditDistance mCheck;

  private byte[] mLastTemplate;
  private int mLastTemplateStart;
  private int mLastTemplateLength;

  private static final int REPEAT = -99;

  private int mCalls = 0;

  private final int[][] mEDCounts;
  private final int[][] mEDScores;

  private static final boolean[][][][] LOWCOMPLEXITY = new boolean[16][16][16][16];

  private final SeedPositions[] mSeeds = new SeedPositions[10000];
  private int mNumSeeds;
  private int[] mWorkspace;

  private final int mGapOpenPenalty;
  private final int mGapExtendPenalty;

  /**
   * Default constructor
   * @param ngsParams supplies other penalty configuration
   * @param nullOnLowComplexity true if this should return null for any alignments which trigger the low complexity filter. This is more accurate, but knocks out lots of alignments.
   */
  SeededAligner(NgsParams ngsParams, boolean nullOnLowComplexity) {
    for (int i = 0; i < mSeeds.length; ++i) {
      mSeeds[i] = new SeedPositions();
    }
    mEDCounts = new int[10][];
    mEDScores = new int[10][];
    for (int i = 0; i < mEDCounts.length; ++i) {
      mEDCounts[i] = new int[10];
      mEDScores[i] = new int[10];
    }
   //resetTemplate();
    mGapOpenPenalty = ngsParams.gapOpenPenalty();
    mGapExtendPenalty = ngsParams.gapExtendPenalty();

    mED = new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false);
    mFixedStart = mED;
    mFixedEnd = mED;
//    mCheck = mED;

    initLowComplexityArray();
    mNullOnLowComplexity = nullOnLowComplexity;
  }

  //TODO: evaluate possible better low complexity values/algorithm
  private void initLowComplexityArray() {
    final int[] count = new int[4];
    for (int i = 0; i <= 15; ++i) {
      for (int j = 0; j <= 15; ++j) {
        for (int k = 0; k <= 15; ++k) {
          for (int l = 0; l <= 15; ++l) {
            count[0] = i;
            count[1] = j;
            count[2] = k;
            count[3] = l;
            for (int m = 0; m < 3; ++m) {
              for (int n = m + 1; n < 4; ++n) {
                if (count[m] < count[n]) {
                  final int tmp = count[m];
                  count[m] = count[n];
                  count[n] = tmp;
                }
              }
            }
            if ((count[0] >= 13) || (count[0] + count[1] >= 12)) {      //TODO 13 and 12 are more or less randomly picked - should come up with reasoned numbers
              LOWCOMPLEXITY[i][j][k][l] = true;
            }
          }
        }
      }
    }

  }

  /**
   * Returns the seeds
   * @return the list of seeds
   */
  /*public SeedPositions[] getSeeds() {
    return mSeeds;
  }*/

  /**
   * Stores a value into a simple hash table
   *
   * @param value the value that is stored
   * @param pos the position in the template
   */
  private void store(final long value, final int pos) {
    //System.err.println("store: "+value+" "+pos);
    int index = (int) (value & mHashtableSize);
    while (mIndex[index] != -1) {
      if (mIndex[index] == value) {
        mPos[index] = REPEAT;
        return;
      }
      ++index;
      if (index >= mIndex.length) {
        index = 0;
      }
    }

    mIndex[index] = value;
    mPos[index] = pos;
  }

  /**
   * Searches for a value in the hash table
   *
   * @param value the hash value to find
   * @return the position (&gt;=0), not found (-1) or duplicate (-99)
   */
  private int lookup(final long value) {
    // System.err.println("find: " + value);
    int index = (int) (value & mHashtableSize);
    while ((mIndex[index] != -1) && (mIndex[index] != value)) {
      ++index;
      if (index >= mIndex.length) {
        index = 0;
      }
    }
    if (mIndex[index] == value) {
      return mPos[index];
    }
    return -1;
  }

  /**
   * Marks the internal template for clearing in the future.
   */
  /*public void resetTemplate() {
    mLastTemplate = null;
    mLastTemplateStart = -1;
    mLastTemplateLength = -1;
    mHashtableSize = -1;
  }*/

  private void allocateWorkspace(int rlen) {
    mHashtableSize = (int) (MathUtils.ceilPowerOf2(rlen * 2L) - 1);
    mIndex = new long[mHashtableSize + 1];
    mPos = new int[mHashtableSize + 1];
    mWorkspace = new int[rlen * 2 + 255]; // why 255?
  }

  /**
   * fast alignment method
   *
   * @param read read sequence
   * @param rlen read length
   * @param template the template
   * @param zeroBasedStart start position
   * @param maxScore maximum edit distance
   * @param maxShift maximum shift for the start or end positions
   * @param cgLeft ignored
   * @return alignment score
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    ++mCalls;
    if (mWorkspace != null) {
      mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
      mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = 0;
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
    }
    if (zeroBasedStart < 0) {
//      if (template.length > 1000) {
//        Diagnostic.developerLog("SeededAligner: -ve start (rejecting): len= " + rlen + " start/end= " + zeroBasedStart + "/" + template.length + " maxScore= " + maxScore);
//      }
      return null; // skipping if really negative start
    } else if (rlen < NTHASHSIZE) {
      return null; // too short to find a hash
    }

    if ((mLastTemplate != template) || (mLastTemplateStart != zeroBasedStart) || (mLastTemplateLength != rlen) || (rlen * 2 > mHashtableSize)) {
      mLastTemplate = template;
      mLastTemplateStart = zeroBasedStart;
      mLastTemplateLength = rlen;
      if (rlen * 2 > mHashtableSize) {
        allocateWorkspace(rlen);
      }
      Arrays.fill(mIndex, -1);
      Arrays.fill(mPos, -1);
      if (!setupTemplateIndex(zeroBasedStart, template, rlen)) {
        return null;
      }
    }
    int state;
    mNumSeeds = 0;
    boolean insameseed = false;
    final int maxOffset = Math.min(maxShift, maxScore);
    // now check the original read
    int fill = 0;
    long rolling = 0;
    for (int i = 0; i < rlen; ++i) {
      final byte b = read[i];
      if (b < 1 || b > 4) {
        fill = 0;
        rolling = 0;
        insameseed = false;
      } else {
        ++fill;
        rolling = (((rolling << 2) + b - 1)) & MAXHASHVALUE;
      }
      // if we have enough bits for lookup, then lookup...
      // if it's not found or over a distance then say we're not in a seed region
      if (fill >= NTHASHSIZE) {
        state = SEEDTYPE;
        final int ret = lookup(rolling);
        if (ret == REPEAT) { // we'll ignore repeated regions in the template
          state = GAPTYPE;
          insameseed = false;
        } else {
          if ((ret == -1) || (Math.abs(ret - (i - NTHASHSIZE + 1 + zeroBasedStart)) > maxOffset)) {
            state = GAPTYPE;
            insameseed = false;
          }
        }
        // if we have found a seed that might be an extension, but it's not in the same diagonal then it's not the same seed
        if ((state == SEEDTYPE) && (mNumSeeds > 0) && insameseed) {
          final int dy = ret - mSeeds[mNumSeeds - 1].mY1;
          final int dx = (i - NTHASHSIZE + 1) - mSeeds[mNumSeeds - 1].mX1;
          if (dx != dy) {
            insameseed = false;
          }
        }
        // if a seed and a new seed and not a repeat then make a seed region
        if (!insameseed && (state == SEEDTYPE) && (ret != REPEAT)) {
          // we only insert a seed
          mSeeds[mNumSeeds].mX1 = i - NTHASHSIZE + 1;
          mSeeds[mNumSeeds].mX2 = mSeeds[mNumSeeds].mX1 + NTHASHSIZE;
          mSeeds[mNumSeeds].mY1 = ret;
          mSeeds[mNumSeeds].mY2 = mSeeds[mNumSeeds].mY1 + NTHASHSIZE;
          mSeeds[mNumSeeds].mType = state;
          insameseed = true;
          ++mNumSeeds;
          if (mNumSeeds >= mSeeds.length) { // far too many seeds
            return null;
          }
        } else {
          if ((state == SEEDTYPE) && (mNumSeeds > 0)) { // extend
            mSeeds[mNumSeeds - 1].mX2++;
            mSeeds[mNumSeeds - 1].mY2++;
          }
        }
      }
    }
    if (mNumSeeds == 0) {
      return null; // check if 0 seeds
    }
//     System.err.println();    SeedPositions.dumpSeeds(mSeeds, mNumSeeds, read, template);
    if (!SeedPositions.seedIntegrity(mSeeds, mNumSeeds, read, rlen, template, zeroBasedStart)) {
      return null;
    }
    removeOverlaps(rlen, template.length);
    if (!SeedPositions.seedIntegrity(mSeeds, mNumSeeds, read, rlen, template, zeroBasedStart)) {
      return null;
    }
    if (!SeedPositions.overlapIntegrity(mSeeds, mNumSeeds, read, rlen, template, zeroBasedStart)) {
      return null;
    }
//    System.err.println();    SeedPositions.dumpSeeds(mSeeds, mNumSeeds, read, template);

    if (mNumSeeds == 0) {
      return null; // check again after cleaning up
    }
    /*
     * TODO possibly only check this if one or both ends aren't tied? A 15 long, low complexity match at both ends of the read and also within
     * maxshift is fairly good evidence of a decent hit... something like if (mSeeds[0].mX1 != 0 || mSeeds[0].mY1 != zeroBasedStart || mSeeds[mNumSeeds - 1].mX2 != rlen || mSeeds[mNumSeeds - 1].mY2 != zeroBasedStart + rlen) {
     * (although should take maxshift into account?
     */
    if (numCovered() < rlen / 2) { // not enough covered, so we have little confidence in the seeds being good
      return null;
    }
    final int score = bestPath(read, rlen, template, zeroBasedStart, maxScore, maxShift);
    if (score == Integer.MAX_VALUE) {
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
      mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
      mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
//    } else if (CHECK_ALIGNMENTS_EVERY_X > 0) {
//      testAlignmentAgreement(read, rlen, template, zeroBasedStart, maxScore, maxShift);
    }
    return mWorkspace;
  }

  /**
   * Sets up an index for the template
   * @param zeroBasedStart the start position
   * @param template the template
   * @param rlen the read length
   * @return false if this alignment should be aborted (failure during index creation)
   */
  private boolean setupTemplateIndex(int zeroBasedStart, byte[] template, int rlen) {
    int fill = 0;
    long rolling = 0;
    final int end = Math.min(zeroBasedStart + rlen, template.length);
    final int[] hist = new int[5];
    for (int i = Math.max(0, zeroBasedStart); i < end; ++i) { // only iterate across the template
      final byte b = template[i];
      if (b < 1 || b > 4) { // if it's an N or something weird reset
        fill = 0;
        rolling = 0;
        hist[0] = 0; hist[1] = 0; hist[2] = 0; hist[3] = 0; hist[4] = 0;
      } else { // otherwise count the bit and stream the hash
        if (i - NTHASHSIZE >= zeroBasedStart) {
          if (hist[template[i - NTHASHSIZE]] > 0) {
            hist[template[i - NTHASHSIZE]]--;
          }
        }
        hist[b]++;
        ++fill;
        rolling = (((rolling << 2) + b - 1)) & MAXHASHVALUE;
      }
      if (fill >= NTHASHSIZE) { // if we have enough bits for hashing, then store
        if (!LOWCOMPLEXITY[hist[1]][hist[2]][hist[3]][hist[4]]) {
          store(rolling, i - NTHASHSIZE + 1); // Who knows what these values should be, but NTHASHSIZE is involved. It's a low complexity filter...
        } else if (mNullOnLowComplexity) {
//           System.err.println("**REJECTING**: " + hist[1] + " " + hist[2] + " " + hist[3] + " " + hist[4]);
//           System.err.println(DnaUtils.bytesToSequenceIncCG(Arrays.copyOfRange(template, i - NTHASHSIZE + 1, i)));
          //this is low complexity. This aligner isn't guaranteed to find correct alignments when hashes are discarded, so terminate.
          return false;
        }
      }
    }
    return true;
  }

  /*
  @JumbleIgnore
  private void testAlignmentAgreement(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift) {
    if ((mCalls & CHECK_ALIGNMENTS_EVERY_X) == 0) {  // every x then try a slow version
      final int[] check = mCheck.calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, false);
      if (check[ActionsHelper.ALIGNMENT_SCORE_INDEX] < mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX]) {
        Diagnostic.developerLog("ED < Seeded Aligner: read=" + DnaUtils.bytesToSequenceIncCG(read)
            + " ed= " + check[ActionsHelper.ALIGNMENT_SCORE_INDEX] + " "
            + DnaUtils.bytesToSequenceIncCG(template, ActionsHelper.zeroBasedTemplateStart(check), Math.min(template.length, rlen + 10))
            + " " + ActionsHelper.toString(check)
            + " seeded= " + mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] + " "
            + DnaUtils.bytesToSequenceIncCG(template, ActionsHelper.zeroBasedTemplateStart(mWorkspace), Math.min(template.length, rlen + 10))
            + " " + ActionsHelper.toString(mWorkspace));
      }
    }
  }*/

  /**
   * Calculate the best path, using gaps and edit distance when required.
   *
   * @param read the original read
   * @param rlen read length
   * @param template template sequence
   * @param zeroBasedStart start in the template
   * @param maxScore maximum score for early termination
   * @param maxShift max shift allowed for start/end positions
   * @return best path
   */
  private int bestPath(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift) {
    int x1, y1;
    int x2, y2;

    if (mNumSeeds == 0) {
      throw new RuntimeException("should be > 0 seeds");
    }

    mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
    mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = 0;
    mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart + rlen;

    int curScore = 0;

    // then the last seed
    if (mNumSeeds > 0) {
      if ((mSeeds[mNumSeeds - 1].mX2 < rlen) || (mSeeds[mNumSeeds - 1].mY2 < template.length)) {
        x1 = mSeeds[mNumSeeds - 1].mX2;
        if (x1 != rlen) {
          // then the gap to the start
          final int[] res = mFixedStart.calculateEditDistanceFixedStart(read, x1, rlen, template, mSeeds[mNumSeeds - 1].mY2, maxScore, maxShift);
          ActionsHelper.prepend(mWorkspace, res);
          mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = mSeeds[mNumSeeds - 1].mY1;

        }
      }
      curScore = mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX];
    }

    // then the middle seeds
    for (int i = mNumSeeds - 1; i >= 1; --i) {
      // first add the seed
      ActionsHelper.prepend(mWorkspace, mSeeds[i].mX2 - mSeeds[i].mX1, ActionsHelper.SAME, 0);
      // then add the gap
      x1 = mSeeds[i - 1].mX2;
      x2 = mSeeds[i].mX1;
      y1 = mSeeds[i - 1].mY2;
      y2 = mSeeds[i].mY1;

      updateEDCounts((x2 - x1) / 10, (y2 - y1) / 10);

      final int xd = x2 - x1;
      final int yd = y2 - y1;
      final long currentScore = curScore; // to avoid integer overflow

      if ((xd != 0) || (yd != 0)) {
        if ((xd == 0) && (yd != 0)) {
          if (currentScore + (yd * (long) mGapExtendPenalty) > maxScore) {
            return Integer.MAX_VALUE;
          }
          ActionsHelper.prepend(mWorkspace, yd, ActionsHelper.DELETION_FROM_REFERENCE, mGapOpenPenalty + yd * mGapExtendPenalty);
        } else if ((yd == 0) && (xd != 0)) {
          if (currentScore + (xd * (long) mGapExtendPenalty) > maxScore) {
            return Integer.MAX_VALUE;
          }
          ActionsHelper.prepend(mWorkspace, xd, ActionsHelper.INSERTION_INTO_REFERENCE, mGapOpenPenalty + xd * mGapExtendPenalty);
        } else {
          final int[] res = mED.calculateEditDistanceFixedBoth(read, x1, x2, template, y1, y2, maxScore - curScore, maxShift);
          if (currentScore + res[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxScore) {    //check this here - if it's integer.maxvalue we don't want to add to it in the next step.
            return Integer.MAX_VALUE;
          }
          ActionsHelper.prepend(mWorkspace, res);
          curScore = mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX];
          if ((curScore >= 0) && (curScore <= maxScore)) {
            updateEDScores((x2 - x1) / 10, (y2 - y1) / 10, curScore);
          }
        }
      }
      if (curScore > maxScore) {
        return Integer.MAX_VALUE;
      }
    }

    // up to the first seed
    mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = mSeeds[0].mY1;

    if (mNumSeeds == 0) {
      x2 = rlen;
      y2 = zeroBasedStart + rlen;
    } else {
      x2 = mSeeds[0].mX1;
      y2 = mSeeds[0].mY1;
    }

    if (mNumSeeds > 0) {
      // first add the seed
      ActionsHelper.prepend(mWorkspace, mSeeds[0].mX2 - mSeeds[0].mX1, ActionsHelper.SAME, 0);
      mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = mSeeds[0].mY1;

      if (mSeeds[0].mX1 != 0) {
        //offset of the first seed between read and template
        final int offset = y2 - zeroBasedStart - x2;
        // then the gap to the start
        final int[] res = mFixedEnd.calculateEditDistanceFixedEnd(read, 0, x2, template, offset + zeroBasedStart, y2, maxScore - curScore, maxShift + Math.abs(offset));
        ActionsHelper.prepend(mWorkspace, res);
      }
    }

    curScore = mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX];

    //maxscore exceeded OR if difference in start positions > maxShift, return maxint.
    if (curScore > maxScore || Math.abs(zeroBasedStart - mWorkspace[ActionsHelper.TEMPLATE_START_INDEX]) > maxShift) {
      return Integer.MAX_VALUE;
    }

//    final ActionsValidator av = new ActionsValidator(1);
//    av.setVerbose(true);
//    if (!av.isValid(mWorkspace, read, rlen, template)) {
//      Diagnostic.developerLog(" read:  " + DnaUtils.bytesToSequenceIncCG(read, 0, rlen));
//      Diagnostic.developerLog(" tmpl:  " + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart, rlen));
//      Diagnostic.developerLog(" sd #:  " + mNumSeeds);
//      Diagnostic.developerLog(" seeds: " + SeedPositions.dumpSeedsAsString(mSeeds));
//      Diagnostic.developerLog(" Invalid action: " + ActionsHelper.toString(mWorkspace));
//      throw new RuntimeException("Invalid action reconstruction");
//    }
//
    return curScore;
  }

  private void updateEDCounts(final int xd, final int yd) {
    if ((xd < 10) && (yd < 10)) {
      mEDCounts[xd][yd]++;
    }
  }

  private void updateEDScores(final int xd, final int yd, int score) {
    if ((xd < 10) && (yd < 10)) {
      mEDScores[xd][yd] += score;
    }
  }

  private int numCovered() {
    int num = 0;
    for (int k = 0; k < mNumSeeds; ++k) {
      num += mSeeds[k].mX2 - mSeeds[k].mX1;
    }
    return num;
  }

  private void deleteSeed(final int i) {
    for (int k = i; k < mNumSeeds - 1; ++k) {
      mSeeds[k].mX1 = mSeeds[k + 1].mX1;
      mSeeds[k].mX2 = mSeeds[k + 1].mX2;
      mSeeds[k].mY1 = mSeeds[k + 1].mY1;
      mSeeds[k].mY2 = mSeeds[k + 1].mY2;
      mSeeds[k].mType = mSeeds[k + 1].mType;
    }
    mSeeds[mNumSeeds].mX1 = -99999999;
    mSeeds[mNumSeeds].mX2 = -99999999;
    mSeeds[mNumSeeds].mY1 = -99999999;
    mSeeds[mNumSeeds].mY2 = -99999999;
    --mNumSeeds;
  }

  /**
   * Remove the overlaps between seeds
   *
   * @param readlen read length
   * @param templatelen template length
   */
  private void removeOverlaps(final int readlen, final int templatelen) {
    // check seeds, what this means is that a seed can't overlap in the x/y direction with the previous seed region
    // truncate seeds which overlap within a small offset of the diagonal, or if they overlap perfectly
    // then delete the big crazy overlaps

    for (int i = 0; i < mNumSeeds - 1; ++i) {
      if ((mSeeds[i].mX2 > mSeeds[i + 1].mX1) || (mSeeds[i].mY2 > mSeeds[i + 1].mY1)) {
//        System.err.println("iY " + mSeeds[i].mY1 + "->" + mSeeds[i].mY2);
//        System.err.println("iX " + mSeeds[i].mX1 + "->" + mSeeds[i].mX2);
//        System.err.println("i+1Y " + mSeeds[i + 1].mY1 + "->" + mSeeds[i + 1].mY2);
//        System.err.println("i+1X " + mSeeds[i + 1].mX1 + "->" + mSeeds[i + 1].mX2);
        if (mSeeds[i + 1].mY1 <= mSeeds[i].mY1 && mSeeds[i + 1].mY2 < mSeeds[i].mY2) {
          //if truncating the seeds would result in a seed with -ve length, just ignore it.
          //Both this and the preceeding seed should be deleted in the next loop.
          continue;
        } else {
          final int dx = Math.abs(mSeeds[i].mX2 - mSeeds[i + 1].mX1);
          final int dy = Math.abs(mSeeds[i].mY2 - mSeeds[i + 1].mY1);

          if (((dx <= 4) && (dy <= 4)) || (dx == dy)) {
  //        System.err.println("truncate a little?" + mSeeds[i] + " " + mSeeds[i + 1] + " " + dx + " " + dy);
            final int move = Math.max(dx, dy);
            mSeeds[i + 1].mX1 += move;
            mSeeds[i + 1].mY1 += move;
          }
        }
      }
    }

    boolean deleted;
    do {
      deleted = false;
      for (int i = 0; i < mNumSeeds - 1; ++i) {
        if ((mSeeds[i].mX2 > mSeeds[i + 1].mX1) || (mSeeds[i].mY2 > mSeeds[i + 1].mY1)) {

          // delete both
          deleteSeed(i + 1);
          deleteSeed(i);
          --i;
          deleted = true;
        }
      }
    } while (deleted);
    // then clean up little overhangs
    for (int i = 0; i < mNumSeeds; ++i) {
      if ((mSeeds[i].mX1 == mSeeds[i].mX2) && (mSeeds[i].mY1 == mSeeds[i].mY2)) {
        // null seed
        deleteSeed(i);
        --i;
      }
    }
  }



  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateStart, int templateEnd, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(final byte[] read, final int readStartPos, final int readEndPos, final byte[] template,
      final int templateExpectedStartPos, final int templateEndPos, final int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  @JumbleIgnore
  public void logStats() {
    final StringBuilder sb = new StringBuilder();
    sb.append("SeededAligner stats").append(StringUtils.LS);
    sb.append("  calls: ").append(mCalls).append(StringUtils.LS);
    for (int i = 0; i < 10; ++i) {
      final int realval = i * 10;
      final String correction = String.format("%6d", realval);
      sb.append(correction).append(" ");
    }
    sb.append(StringUtils.LS); // for the table
    for (int i = 0; i < 10; ++i) {
      sb.append("------");
    }
    sb.append(StringUtils.LS); // for the table
    for (int j = 0; j < 10; ++j) {
      for (int i = 0; i < 10; ++i) {
        final String correction = String.format("%6d", mEDCounts[j][i]);
        sb.append(correction).append(" ");
      }
      sb.append(StringUtils.LS);
    }
    sb.append(StringUtils.LS); // for the table
    for (int j = 0; j < 10; ++j) {
      for (int i = 0; i < 10; ++i) {
        final int num =  mEDCounts[j][i] > 0 ? mEDScores[j][i] / mEDCounts[j][i] : 0;
        final String correction = String.format("%6d", num);
        sb.append(correction).append(" ");
      }
      sb.append(StringUtils.LS);
    }
    Diagnostic.developerLog(sb.toString()); // log in one call
  }
}
