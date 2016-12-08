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
package com.rtg.ngs.blocking;

import java.io.Closeable;
import java.util.Arrays;

import com.rtg.util.License;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Maintain a count of number of records at the best score and the
 * next best score for each read. Use these to calculate an
 * approximate mapping quality score using the same calculation that
 * <code>BWA</code> uses.
 *
 */
public class MapQScoringReadBlocker implements Closeable {

  protected static final int MAX_COUNT = Short.MAX_VALUE * 2 + 1;

  protected static final int MAX_SCORE = Byte.MAX_VALUE * 2 + 1;

  protected static final int NEXT_LIMIT = 255;

  private static final int[] G_LOG_N;
  static {
    G_LOG_N = new int[NEXT_LIMIT + 1];
    for (int i = 1; i < G_LOG_N.length; ++i) {
      G_LOG_N[i] = (int) (4.343 * Math.log(i) + 0.5);
      //System.err.println("G_LOG_N[" + i + "]=" + G_LOG_N[i] + " 23-=" + ((23 < G_LOG_N[i]) ? 0 : 23 - G_LOG_N[i]));
    }
  }

  /**
   * Based on BWA approximate MAPQ score calculation
   *
   * @param bestCount number of hits at the best score
   * @param nextCount number of hits at the second best score
   * @param bestScore alignment score of the best score
   * @param maxScore user specified maximum alignment score
   * @return a MAPQ score
   */
  public static int approxMapQ(final int bestCount, final int nextCount, final int bestScore, final int maxScore) {

    // Score for unmapped reads (probability that the read does not come from the reference)
    if (bestCount == 0) {
      return 23;
    }

    // Since error prob must be greater than 1/bestCount in this case
    if (bestCount > 1) {
      return 0;
    }

    // Maybe in this case there will be no second best count (they
    // would be worse than maxAlignmentScore), so downrate??
    if (bestScore == maxScore) {
      return 25;
    }

    // Probability of unique hit being correct when no second best hits?
    if (nextCount == 0) {
      return 37;
    }

    // Otherwise decrease quality according to number of next best hits
    final int n = (nextCount >= NEXT_LIMIT) ? NEXT_LIMIT : nextCount;
    return (23 < G_LOG_N[n]) ? 0 : 23 - G_LOG_N[n];
  }

  private final String mTitle;
  private final int mThreshold;  // Count limit for best score
  private final byte[] mScore1; // Current best score treated as unsigned here
  private final byte[] mScore2; // Current next best score treated as unsigned here
  private final short[] mCount1; // Count at current score treated as unsigned here
  private final short[] mCount2; // Count at next best score treated as unsigned here

  /**
  * Copy constructor for making a non-synchronized version from a synchronized one
   * @param source the source MapQScoringReadBlocker
   */
  public MapQScoringReadBlocker(MapQScoringReadBlocker source) {
    mTitle = source.mTitle;
    mThreshold = source.mThreshold;
    mScore1 = source.mScore1;
    mScore2 = source.mScore2;
    mCount1 = source.mCount1;
    mCount2 = source.mCount2;
  }

  /**
   * Creates a counter for <code>count</code> records blocking at <code>threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold score in range 1 to 255
   * @param title a title to use during logging
   */
  public MapQScoringReadBlocker(final int count, final int threshold, final String title) {
    // we add one to the threshold because the way this scoring blocker is
    // used means that it discards reads with count >= the threshold.
    if (threshold > MAX_COUNT || threshold < 1) {
      throw new IllegalArgumentException("threshold must be > 0, not " + threshold);
    }
    mTitle = title;
    mThreshold = threshold + 1;
    mCount1 = new short[count];
    mCount2 = new short[count];
    mScore1 = new byte[count];
    mScore2 = new byte[count];
    Arrays.fill(mScore1, (byte) 0xFF);
    Arrays.fill(mScore2, (byte) 0xFF);
  }

  /**
   * Creates a counter for <code>count</code> records blocking at <code>
   * threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold in range 1 to 255.
   */
  public MapQScoringReadBlocker(final int count, final int threshold) {
    this(count, threshold, "blocked pairings");
  }

  /**
   * The number of records with this read identifier at the best score
   * If this is equal to or greater than the threshold, then it means
   * that a large/infinite number of records have been written.
   * @param r Read identifier
   * @return the number of records written, or the threshold.
   */
  public final int getCount1(int r) {
    return mCount1[r] & MAX_COUNT;
  }

  /**
   * The best score.
   * @param r Read identifier
   * @return the score.
   */
  public final int getScore1(int r) {
    return mScore1[r] & 0xFF;
  }

  /**
   * The number of records with this read identifier at the second
   * best score. If this is equal to or greater than the threshold,
   * then it means that a large/infinite number of records have been
   * written.
   * @param r Read identifier
   * @return the number of records written, or the threshold.
   */
  final int getCount2(int r) {
    return mCount2[r] & MAX_COUNT;
  }

  /**
   * The second best score.
   * @param r Read identifier
   * @return the score.
   */
  public final int getScore2(int r) {
    return mScore2[r] & 0xFF;
  }

  /**
   * Return the score which can be used for the early termination of alignments
   * @param r Read identifier
   * @return the score
   */
  public int getTerminationScore(int r) {
    if (getCount1(r) > 1) {
      return getScore1(r);
    }
    return getScore2(r);
  }

  /**
   * Increment the count for the given read.
   *
   * @param r read number
   * @param score the score to associate with the read
   * @return the count at the current best score or -1 if we did not increment at the best score
   */
  public int increment(final int r, int score) {
    //System.err.println("increment(" + r + "," + score + ")");
    final int currentScore = mScore1[r] & 0xFF;
    final int currentScore2 = mScore2[r] & 0xFF;
    if (score < currentScore) {       // New best score, old best becomes second best
      mScore2[r] = mScore1[r];
      mCount2[r] = mCount1[r];
      mScore1[r] = (byte) ((score > MAX_SCORE) ? 0xFF : (score & 0xFF));
      mCount1[r] = 1;
      return 1;
    } else if (score == currentScore) { // Another hit at the best score
      if (mCount1[r] != -1) {
        mCount1[r]++;
      }
      return getCount1(r);
    } else if (score < currentScore2) { // New second best score
      mScore2[r] = (byte) ((score > MAX_SCORE) ? 0xFF : (score & 0xFF));
      mCount2[r] = 1;
    } else if (score == currentScore2) { // Another hit at second best score
      if (mCount2[r] != -1) {
        mCount2[r]++;
      }
    } // score > currentScore2 is a no-op
    return -1;
  }

  /**
   * Returns a MAPQ value for a mated read.
   *
   * @param r read id
   * @return a MAPQ value
   */
  public final int getMatedMapQ(final int r) {
    return Math.min(getMapQ(r) * 3 / 2, 254);   // For mated assume roughly 1.5 quality, but don't exceed max allowed
  }

  /**
   * Returns a MAPQ value for an unmated read.
   *
   * @param r read id
   * @return a MAPQ value
   */
  public final int getMapQ(final int r) {
    final int bestCount = mCount1[r] & MAX_COUNT;
    final int nextCount = mCount2[r] & MAX_COUNT;
    final int bestScore = mScore1[r] & 0xFF;
    return approxMapQ(bestCount, nextCount, bestScore, MAX_SCORE); // MAX_SCORE should actually be the per read max score derived from the user -e/-E max (incorporating current read length where appropriate)
  }

  /**
   * Check if the specified read is blocked at the best score.
   *
   * @param r read to check
   * @param score the score to check
   * @return true if blocked
   */
  public boolean isBlocked1(final int r, int score) {
    final int currentScore = mScore1[r] & 0xFF;
    return (score > currentScore)
    || ((score == currentScore) && (mCount1[r] & MAX_COUNT) >= mThreshold);
  }

  /**
   * Given read has too many hits
   * @param r read id
   * @return true if read is blocked
   */
  public boolean isBlocked(final int r) {
    return (mCount1[r] & MAX_COUNT) >= mThreshold;
  }

  /**
   * Check if the specified read is blocked at the second best score.
   *
   * @param r read to check
   * @param score the score to check
   * @return true if blocked
   */
  public boolean isBlocked2(final int r, int score) {
    final int currentScore = mScore2[r] & 0xFF;
    return (score > currentScore)
    || ((score == currentScore) && (mCount2[r] & MAX_COUNT) >= NEXT_LIMIT);
  }

  @Override
  public void close() {
    if (License.isDeveloper()) {
      final int[] h = new int[MAX_COUNT + 1];
      for (final int k : mCount1) {
        h[k & MAX_COUNT]++;
      }
      Diagnostic.developerLog("Statistics of " + mTitle);
      long sum = 0;
      for (int k = 0; k < MAX_COUNT + 1; ++k) {
        if (h[k] > 0) {
          sum += h[k];
          final String c = k == MAX_COUNT ? ">= " + MAX_COUNT : String.valueOf(k);
          Diagnostic.developerLog(h[k] + " reads had count " + c);
        }
      }
      Diagnostic.developerLog("Total reads " + sum);
    }
  }

  @Override
  public String toString() {
    return "ScoringReadBlocker";
  }

}
