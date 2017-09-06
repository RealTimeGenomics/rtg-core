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

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Fast aligner for a substitution, an <code>N</code> difference, or a single
 * insert/delete every 10 places.  This is intended for long reads, 100 or greater.
 *
 */
class HopStepEditDistanceLong implements UnidirectionalEditDistance, Integrity {

  /** Minimum number of equal bases we have to find after a MNP/insert/delete. */
  private static final int RESYNC_LENGTH = 7;

  /** Maximum length of each insert/delete/MNP. */
  private static final int MAX_INDEL = 4;

  /** Stop the hop-skip algorithm when we get closer than this to the start. */
  private static final int TOO_NEAR_TO_START = RESYNC_LENGTH + 2 * MAX_INDEL + 1;

  /** Minimum distance between indels */
  private static final int INDEL_SPACING = 21;

  /** A <code>SeedPositions.mType</code> code that means align this region later. */
  private static final int DO_IT_LATER = 99;

  private int mReadLengthCutoff = 50;

  private final UnidirectionalEditDistance mEd;

  /** Gap opening penalty. */
  private final int mGapOpen;
  private final int mGapExtend;
  private final int mSubstitionPenalty;
  private final int mUnknownsPenalty;

  private int mReadLen;
  private byte[] mRead;
  private byte[] mTemplate;
  private int mNumSeeds;

  private int[] mWorkspace;

  private SeedPositions[] mSeeds;

  /** A spare one, for <code>findSeed</code> to record its results */
  private SeedPositions mEqSeed = new SeedPositions();

  private int[] mHistogram = new int[10];
  private static final int HISTOGRAM_BIN = 10;
  private int mStatsTotal;

  /** These counters, plus the histogram, should add up to total. */
  private int mStatsSuccess, mStatsOffTemplate, mStatsCannotSync, mStatsMaxShift, mStatsIndelsTooClose,
    mStatsEndTooHard, mStatsMiddleTooHard, mStatsMiddleMerged, mStatsStartTooHard, mStatsStartMoved;

  /** These count the number of times we call the full aligner. */
  private int mStatsEndIsComplex, mStatsMiddleIsComplex, mStatsStartIsComplex;

  /**
   * New fast aligner for highly similar sequences.
   * @param params {@link NgsParams} for current run
   */
  protected HopStepEditDistanceLong(NgsParams params) {
    mSeeds = new SeedPositions[10];
    for (int i = 0; i < mSeeds.length; ++i) {
      mSeeds[i] = new SeedPositions();
    }
    mNumSeeds = 0;
    mGapOpen = params.gapOpenPenalty();
    mGapExtend = params.gapExtendPenalty();
    mSubstitionPenalty = params.substitutionPenalty();
    mUnknownsPenalty = params.unknownsPenalty();
    mEd = new GotohEditDistance(params);
    integrity();
  }

  void setReadLengthCutoff(int cutoff) {
    mReadLengthCutoff = cutoff;
  }

  private void checkWorkspace(final int len) {
    if (mWorkspace == null || mWorkspace.length < len) {
      mWorkspace = new int[len];
    }
  }

  private void incrHistogram(final int x) {
    final int pos = x / HISTOGRAM_BIN;
    if (pos >= mHistogram.length) {
      mHistogram = Arrays.copyOf(mHistogram, pos + 10);
    }
    mHistogram[pos]++;
  }

  /**
   * Search around in the template to find a perfect match against the read.
   * This tries <code>mTries</code> different positions along the read,
   * moving down <code>mJumpSize</code> each time.
   * For each read position, it tries to match it against the corresponding
   * segment of the template, shifting the template position around
   * by up to +/- <code>maxY - y</code>.
   *
   * Precondition: the comparison region of the read is completely inside the read.
   *
   * @param x end of the read segment to match
   * @param y approximate corresponding position in the template
   * @param maxX rightmost position in read that we are allowed to match (inclusive)
   * @param maxY rightmost position in template that we are allowed to match (inclusive)
   * @param windowLen the number of bases that must match
   * @param tries number of times to jump down along the read trying to find a seed
   * @param jumpSize the size of each jump
   * @return the seed that was found, else null
   */
  SeedPositions findSeed(int x, int y, int maxX, int maxY, int windowLen, int tries, int jumpSize) {
    assert windowLen - 1 <= x && x < mReadLen;
    assert x <= maxX;
    assert y < maxY;
    final int maxShift = maxY - y + 2; // two extra, to make really sure the match is unique
    // move down the read, looking for a perfect equality region.
    for (int j = 0; j < tries; ++j) {
      final int xEnd = x - j * jumpSize;
      final int yEnd = y - j * jumpSize;
      if (xEnd + 1 - windowLen < 0) {
        continue;
      }
      // find best end position.  Delta = 0, +1, -1, +2, -2, +3, -3, +4,...
      int matchCount = 0; // the number of ways we can resynchronize
      int bestDelta = 0; // the delta of the first good match
      for (int delta = 0; delta <= maxShift; delta = -delta + (delta <= 0 ? 1 : 0)) {
        final int ydelta = yEnd + delta;
        if (ydelta >= mTemplate.length || ydelta - (windowLen - 1) < 0) {
          continue;
        }
        boolean matches = true;
        for (int i = 0; i < windowLen; ++i) {
          final byte r = mRead[xEnd - i];
          final byte t = mTemplate[ydelta - i];
          if (r * t == 0) {
            return null;  //just say we can't do Ns. I can't figure out how to hack this algorithm to support unknowns with different penalties to mismatches.
          }
          if (r != t) {
            matches = false;
            break; // since we require perfect matches
          }
        }
        if (matches) {
          ++matchCount;
          if (matchCount == 1) {
            // remember this match
            bestDelta = delta;
          } else if (matchCount > 1) {
            break;  // it is ambiguous how to resynchronize here
          }
        }
      }
      if (matchCount == 1 && Math.abs(bestDelta) <= maxShift - 2) {
        // there is a unique way of resynchronizing
        mEqSeed.mType = ActionsHelper.SAME;
        mEqSeed.mX1 = xEnd + 1 - windowLen;
        mEqSeed.mX2 = xEnd + 1;
        mEqSeed.mY1 = yEnd + bestDelta + 1 - windowLen;
        mEqSeed.mY2 = yEnd + bestDelta + 1;
//        System.err.println("extending seed " + mEqSeed.toString());
        extendSeedRight(mEqSeed, maxX + 1, Math.min(maxY + 1, mTemplate.length));
//        System.err.println("extended seed " + mEqSeed.toString());
        return mEqSeed;
      }
    }
    return null;
  }

  /**
   * Extends the given equality region as far as it is equal,
   * but limited by the two given bounds.
   *
   * @param seed the existing seed to extend
   * @param maxX the maximum bound on the read position (exclusive)
   * @param maxY maximum bound on the template position (exclusive)
   */
  void extendSeedRight(final SeedPositions seed, final int maxX, final int maxY) {
    assert maxX <= mReadLen;
    assert maxY <= mTemplate.length;
    int xEnd = seed.mX2;
    int yEnd = seed.mY2;
    while (xEnd < maxX && yEnd < maxY && !different(mRead[xEnd], mTemplate[yEnd])) {
      ++xEnd;
      ++yEnd;
    }
    seed.mX2 = xEnd;
    seed.mY2 = yEnd;
  }

  private int[] failure() {
    mRead = null;
    mTemplate = null;
    mNumSeeds = 0;
    return null;
  }

  private void initEd(byte[] read, int rlen, byte[] template) {
    mReadLen = rlen;
    mRead = read;
    mTemplate = template;
    mNumSeeds = 0;
    ++mStatsTotal;
    mWorkspace[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
    mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] = 0;
  }

  /**
   * Hopping, skipping edit distance.
   *
   * @param read the encoded read
   * @param rlen length of read (assumed that given read is at least this long)
   * @param template the encoded template
   * @param zeroBasedStart start position in the template
   * @param maxScore the maximum alignment score, for early termination
   * @param maxShift the maximum shift allowed for start or end position
   * @param cgLeft ignored
   * @return the actions array, or null if not possible
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    if (rlen < mReadLengthCutoff) { //should only use this aligner for longish reads
      return null;
    }
    checkWorkspace(rlen + ActionsHelper.ACTIONS_START_INDEX + 1);
    initEd(read, rlen, template);
    int x = rlen - 1;
    int y = zeroBasedStart + rlen - 1;
    if (y >= template.length + MAX_INDEL || y < rlen - MAX_INDEL) {
      ++mStatsOffTemplate;
      return failure();
    }
    final int windowSize = Math.min(10, rlen);
    final int tries = 4;
    final int jumpSize = 5;
    final SeedPositions eq = findSeed(x, y, x, y + maxShift, windowSize, tries, jumpSize);
    if (eq == null) {
      ++mStatsCannotSync;
      return failure();
    }
    // use eq as our first equality seed.  This may be a bit before the end of the read.
    mEqSeed = mSeeds[0]; // recycle it for findSeed.
    mSeeds[0] = eq;
    mNumSeeds = 1;
    x = eq.mX1;
    y = eq.mY1;
    mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] = y; // we override this later
    int lastIndelSeed = Integer.MIN_VALUE;
    int windowLength2 = RESYNC_LENGTH;
    int tries2 = 3;
    int[] alignStart = null; // the alignment for the first few nucleotides
    while (x > 0) {
      if (Math.abs(x - (y - zeroBasedStart)) > maxShift) {
        ++mStatsMaxShift;
        return failure();
      }
      boolean allEqual = true;
      while (x > 0) { // this is the inner loop that deals with long strings of equalities
        if (y <= 0) {
          ++mStatsOffTemplate;
          return failure();
        }
        final byte r = read[--x];
        final byte t = template[--y];
        if (r * t == 0) {
          return failure(); //can't handle Ns anymore
        } else if (r != t) {
          allEqual = false;
          break;
        }
      }
      if (allEqual) {
        assert x == 0;
        mSeeds[mNumSeeds - 1].mX1 = x;
        mSeeds[mNumSeeds - 1].mY1 = y;
        break;
      }
      // update the start of the current equality region
      mSeeds[mNumSeeds - 1].mX1 = x + 1;
      mSeeds[mNumSeeds - 1].mY1 = y + 1;
      if (x < TOO_NEAR_TO_START) {      // too near the start, so use a proper aligner to analyze it.
        ++mStatsStartIsComplex;
        alignStart = mEd.calculateEditDistanceFixedEnd(read, 0, x + 1, template, y - x, y + 1, maxScore, maxShift);
        if (alignStart == null || alignStart[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxScore) {
          ++mStatsStartTooHard;
          return failure();
        }
        y = alignStart[ActionsHelper.TEMPLATE_START_INDEX];
        alignStart = ActionsHelper.copy(alignStart);
        break;
      } else if (x < TOO_NEAR_TO_START * 2) {
        // be more careful near the beginning of the read
        windowLength2 = RESYNC_LENGTH * 3 / 2;
        tries2 = tries;
      }
      assert different(read[x], template[y]);
      assert !different(read[x + 1], template[y + 1]);
      final SeedPositions resync = findSeed(x - MAX_INDEL, y - MAX_INDEL, x, y, windowLength2, tries2, RESYNC_LENGTH);
      if (resync == null) { // could not resynchronize
        incrHistogram(x);
        return failure();
      }
      final int eqX = resync.mX2 - 1;
      final int eqY = resync.mY2 - 1;
      assert !different(read[eqX], template[eqY]);
      assert eqX <= x;
      assert eqY <= y;
      final int type;
      final boolean isIndel;
      if (eqX == x) {
        type = ActionsHelper.DELETION_FROM_REFERENCE;
        isIndel = true;
      } else if (eqY == y) {
        type = ActionsHelper.INSERTION_INTO_REFERENCE;
        isIndel = true;
      } else if (x - eqX == y - eqY && x - eqX <= 2) {
          type = ActionsHelper.MISMATCH;
          isIndel = false;
      } else {
        type = DO_IT_LATER; // use a proper aligner later
        isIndel = true; // big MNPs might hide indels, so we merge them with nearby indels.
      }
      if (isIndel && lastIndelSeed == mNumSeeds - 2 && mSeeds[mNumSeeds - 1].xWidth() < INDEL_SPACING) {
        if (mSeeds[lastIndelSeed].xWidth() > INDEL_SPACING) {
          ++mStatsIndelsTooClose;
          return failure(); // that region has probably been merged once already
        }
        ++mStatsMiddleMerged; // we merge these two regions
        mSeeds[lastIndelSeed].mType = DO_IT_LATER;
        mSeeds[lastIndelSeed].mX1 = eqX + 1;
        mSeeds[lastIndelSeed].mY1 = eqY + 1;
        mNumSeeds = lastIndelSeed + 1;
      } else {
        if (isIndel) {
          lastIndelSeed = mNumSeeds; // the one we are about to add.
        }
        addCigar(eqX + 1, eqY + 1, type);
      }
      // check that findSeed has set up the correct equality region
      assert resync.mType == ActionsHelper.SAME;
      assert resync.mX1 <= eqX + 1 - RESYNC_LENGTH;
      assert resync.mX2 == eqX + 1;
      assert resync.mY1 <= eqY + 1 - RESYNC_LENGTH;
      assert resync.mY2 == eqY + 1;
      mEqSeed = mSeeds[mNumSeeds]; // recycle for findSeed
      mSeeds[mNumSeeds++] = resync;
      x = resync.mX1;
      y = resync.mY1;
    } // end of main alignment while loop

    final int[] result = buildAlignment(maxScore, maxShift);
    if (result == null) {
      return failure();
    } else if (alignStart != null) {
      ActionsHelper.prepend(result, alignStart);
    }
    if (Math.abs(zeroBasedStart - y) > maxShift) {
      ++mStatsStartMoved;
      return failure(); // there might be a worse alignment that does not shift so far...
    }
    result[ActionsHelper.TEMPLATE_START_INDEX] = y;
    if (result[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxScore) {       // we are sure the alignment is too bad
//      return failure();   //I don't trust this edit distance, and think this should return null to be safe; but not doing so due to possible performance implications
      result[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
      result[ActionsHelper.ACTIONS_LENGTH_INDEX] = 0;
    }
    ++mStatsSuccess;
    return result;
  }

  private int[] buildAlignment(int maxScore, int maxShift) {
    //System.err.println();    SeedPositions.dumpSeeds(mSeeds, mNumSeeds, mRead, mTemplate);
    if (mSeeds[0].mX2 != mReadLen) {
      ++mStatsEndIsComplex; // We have a SNP/indel at the end.
      final int xEnd = mSeeds[0].mX2;
      final int yEnd = mSeeds[0].mY2;
      //System.err.println("handling end case: " + xEnd + " .. " + rlen);
      final int[] alignEnd = mEd.calculateEditDistanceFixedStart(mRead, xEnd, mReadLen, mTemplate, yEnd, maxScore, maxShift);
      if (alignEnd == null || alignEnd[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxScore) {
        ++mStatsEndTooHard;
        return null;
      }
      ActionsHelper.prepend(mWorkspace, alignEnd);
    }
    for (int i = 0; i < mNumSeeds && mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX] <= maxScore; ++i) {
      //System.err.println("seed[" + i + "] = " + mSeeds[i].toString());
      if (mSeeds[i].mType == ActionsHelper.SAME) {
        final int len = mSeeds[i].mX2 - mSeeds[i].mX1;
        ActionsHelper.prepend(mWorkspace, len, ActionsHelper.SAME, 0);
      } else if (mSeeds[i].mType == ActionsHelper.MISMATCH) {
        final int len = mSeeds[i].mX2 - mSeeds[i].mX1;
        ActionsHelper.prepend(mWorkspace, len, ActionsHelper.MISMATCH, len * mSubstitionPenalty);
      } else if (mSeeds[i].mType == ActionsHelper.INSERTION_INTO_REFERENCE) {
        final int len = mSeeds[i].mX2 - mSeeds[i].mX1;
        ActionsHelper.prepend(mWorkspace, len, ActionsHelper.INSERTION_INTO_REFERENCE, len * mGapExtend + mGapOpen);
      } else if (mSeeds[i].mType == ActionsHelper.DELETION_FROM_REFERENCE) {
        final int len = mSeeds[i].mY2 - mSeeds[i].mY1;
        ActionsHelper.prepend(mWorkspace, len, ActionsHelper.DELETION_FROM_REFERENCE, len * mGapExtend + mGapOpen);
      } else if (mSeeds[i].mType == DO_IT_LATER) {
        ++mStatsMiddleIsComplex;
        final SeedPositions r = mSeeds[i];
//        System.err.println("handling complex region: " + r + " width=" + r.xWidth());
        final int maxMidScore = maxScore - ActionsHelper.alignmentScore(mWorkspace);
        final int[] align = mEd.calculateEditDistanceFixedBoth(mRead, r.mX1, r.mX2, mTemplate, r.mY1, r.mY2, maxMidScore, maxShift);
//        System.err.println(ActionsHelper.toString(align));
        if (align == null || align[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxMidScore) {
          ++mStatsMiddleTooHard;
          return null;
        }
        ActionsHelper.prepend(mWorkspace, align);
      } else {
        Diagnostic.developerLog("HopStep: > 4 type");
        return null;
      }
    }
    return mWorkspace;
  }

  void validateResult(final int zeroBasedStart, final int maxScore) {
    final ActionsValidator av = new ActionsValidator(mGapOpen, mGapExtend, mSubstitionPenalty, mUnknownsPenalty);
    if (!av.isValid(mWorkspace, mRead, mReadLen, mTemplate, maxScore)) {
      //System.err.println(ActionsHelper.toString(mWorkspace));
      final String error = "HopStepEditDistance validation problem: "
        + ActionsHelper.matchCount(mWorkspace) + StringUtils.LS
        + av.getErrorDetails(mWorkspace, mRead, mReadLen, mTemplate, zeroBasedStart);
      Diagnostic.developerLog(error);
      //throw new RuntimeException(error);
    }
  }

  /**
   * True if two bases are not equal, taking n handling into account (n = always different)
   * @param readBase the read value
   * @param templateBase the template value
   * @return true if they are not equal.
   */
  final boolean different(byte readBase, byte templateBase) {
    return readBase != templateBase || readBase == DnaUtils.UNKNOWN_RESIDUE;
  }

  // just for testing
  void setupTestData(final byte[] read, final byte[] template) {
    mRead = read;
    mReadLen = read.length;
    mTemplate = template;
    mNumSeeds = 1;
  }

  private void addCigar(int x, int y, int type) {
    if (mNumSeeds >= mSeeds.length - 2) {
      final int oldLength = mSeeds.length;
      mSeeds = Arrays.copyOf(mSeeds, mNumSeeds * 2);
      for (int i = oldLength; i < mSeeds.length; ++i) {
        mSeeds[i] = new SeedPositions();
      }
    }
    final int prevX = mSeeds[mNumSeeds - 1].mX1;
    final int prevY = mSeeds[mNumSeeds - 1].mY1;
    assert x <= prevX;
    assert y <= prevY;
    assert type != ActionsHelper.DELETION_FROM_REFERENCE || x == prevX;
    assert type != ActionsHelper.INSERTION_INTO_REFERENCE || y == prevY;
    mSeeds[mNumSeeds].mType = type;
    mSeeds[mNumSeeds].mX1 = x;
    mSeeds[mNumSeeds].mX2 = prevX;
    mSeeds[mNumSeeds].mY1 = y;
    mSeeds[mNumSeeds].mY2 = prevY;
    ++mNumSeeds;
  }

  @Override
  public void logStats() {
    final StringBuilder sb = new StringBuilder();
    sb.append(this.toString()).append(" statistics:").append(StringUtils.LS);
    sb.append("Total: ").append(mStatsTotal).append(StringUtils.LS);
    sb.append("Success: ").append(mStatsSuccess).append(StringUtils.LS);
    sb.append("Cannot synch  : ").append(mStatsCannotSync).append(StringUtils.LS);
    sb.append("Off template  : ").append(mStatsOffTemplate).append(StringUtils.LS);
    sb.append("Close Indels  : ").append(mStatsIndelsTooClose).append(StringUtils.LS);
    sb.append("Shifted > maxShift: ").append(mStatsMaxShift).append(StringUtils.LS);
    sb.append("End too hard  : ").append(mStatsEndTooHard).append(" / ").append(mStatsEndIsComplex).append(StringUtils.LS);
    sb.append("Mid too hard  : ").append(mStatsMiddleTooHard).append(" / ").append(mStatsMiddleIsComplex).append(" with ").append(mStatsMiddleMerged).append(" merges").append(StringUtils.LS);
    sb.append("Start too hard: ").append(mStatsStartTooHard).append(" / ").append(mStatsStartIsComplex).append(StringUtils.LS);
    sb.append("Start moved   : ").append(mStatsStartMoved).append(StringUtils.LS);
    for (int i = 0; i < mHistogram.length; ++i) {
      sb.append(i * HISTOGRAM_BIN).append("..: ").append(mHistogram[i]).append(StringUtils.LS);
    }
    Diagnostic.developerLog(sb.toString());
  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(final byte[] read, final int readStartPos, final int readEndPos, final byte[] template, final int templateExpectedStartPos,
      final int templateEndPos, final int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public String toString() {
    // just static information
    return "HopStepEditDistanceLong["
    + "resync_length=" + RESYNC_LENGTH
    + ", max_indel=" + MAX_INDEL
    + ", indel_spacing=" + INDEL_SPACING
    + "]";
  }

  @Override
  public final boolean integrity() {
    Exam.assertTrue(0 <= mNumSeeds && mNumSeeds < mSeeds.length);
    Exam.assertEquals(mRead != null, mNumSeeds > 0);
    if (mNumSeeds > 0) {
      Exam.assertTrue(0 < mReadLen && mReadLen <= mRead.length);
    }
    // check that seeds are contiguous and in reverse order.
    for (int i = 1; i < mNumSeeds; ++i) {
      Exam.assertEquals(mSeeds[i].mX2, mSeeds[i - 1].mX1);
      Exam.assertEquals(mSeeds[i].mY2, mSeeds[i - 1].mY1);
      // and that each seed is well-formed
      Exam.assertTrue(mSeeds[i].mX1 <= mSeeds[i].mX2);
      Exam.assertTrue(mSeeds[i].mX1 <= mSeeds[i].mX2);
      // and non-empty
      Exam.assertTrue((mSeeds[i].mX2 - mSeeds[i].mX1) + (mSeeds[i].mY2 - mSeeds[i].mY1) > 0);
    }
    // Note: the seeds do not necessarily over the whole read, since we may give
    // up aligning part-way through a read, and because the ends are sometimes
    // done using the full aligner instead.
    return false;
  }

  @Override
  public boolean globalIntegrity() {
    return integrity();
  }
}
