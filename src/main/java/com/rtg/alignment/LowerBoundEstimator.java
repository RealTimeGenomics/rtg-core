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

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Simple four nucleotide based lower bound estimator.
 *
 * Gives any read/template with unknowns a score of 0 (unless the read has already been rejected
 * before the unknown regions are encountered.
 *
 */
final class LowerBoundEstimator {

  private final int mWordSize;
  private final int mArraySize;

  private final int[] mOccurs;
  private int[] mWords;
  private final int[] mCount;

  private final int mSubstitutionPenalty;
  private final int mUnknownsPenalty;

  /**
   * Default constructor
   * @param substitutionPenalty the penalty for a substitution
   * @param unknownsPenalty the penalty for an unknown nucleotide
   */
  LowerBoundEstimator(int substitutionPenalty, int unknownsPenalty) {
    this(4, substitutionPenalty, unknownsPenalty);
  }

  /**
   * Default constructor
   * @param wordsize the word size
   * @param substitutionPenalty the penalty for a substitution
   * @param unknownsPenalty the penalty for an unknown nucleotide
   */
  protected LowerBoundEstimator(int wordsize, int substitutionPenalty, int unknownsPenalty) {
    mWordSize = wordsize;
    mUnknownsPenalty = unknownsPenalty;
    mArraySize = (int) Math.pow(4, wordsize) - 1;
    mOccurs = new int[mArraySize + 1];
    mWords = new int[5];
    mCount = new int[wordsize];
    mSubstitutionPenalty = substitutionPenalty;
  }

  private void resizeArrays(final int rlen) {
    if (mWords.length < rlen) {
      mWords = new int[rlen];
    }
  }

  int getCount(final int pos) {
    return mCount[pos];
  }

  /**
   * Lower bound estimator, using a 4-nt lock-step creeping count.
   * @param read the read
   * @param rStartPos where in the read to start (zero based)
   * @param rlen read length
   * @param template the template
   * @param zeroBasedStart start position
   * @param maxScore maximum edit distance
   * @param maxShift the maximum allowed shift in start and end positions
   * @return alignment score
   */
  private int calcLB(byte[] read, int rStartPos, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift) {
    int fill = 0;
    int rolling = 0;

    if (zeroBasedStart >= template.length) {
      Diagnostic.developerLog("calcLB out of range: " + zeroBasedStart + " " + template.length + " " + maxScore);
      return 0; // lower bound of 0
    }

    if (mUnknownsPenalty <= 0) {
      if (zeroBasedStart - maxShift < 0 || zeroBasedStart + rlen + maxShift > template.length) {
        return 0;
      }
    }

    Arrays.fill(mCount, 0);
    Arrays.fill(mOccurs, 0);
    Arrays.fill(mWords, 0);

    resizeArrays(rlen);

    int templatePos = zeroBasedStart < 0 ? 0 : zeroBasedStart;
    // first go maxShift in advance
    final int offset = Math.min(maxShift, maxScore);
    while (templatePos < (zeroBasedStart + offset) && (templatePos - zeroBasedStart) < template.length) {
      byte b = 0;
      //System.err.println("first off: " + templatePos + " fill=" + fill + " z" + zeroBasedStart);
      if ((templatePos >= 0) && (templatePos < template.length)) {
        b = template[templatePos];
        if (b == 0) {
          // if N then return!
          //System.err.println("A");
          return 0;
        }
      }
      //System.err.println("b=" + Utils.getBase(b) + " TP" + templatePos + " " + template.length + " t[0]" + template[0]);
      if (b < 1 || b > 4) { // if it's an N or something weird reset
        return 0;
      } else { // otherwise count the bit and stream the hash
        ++fill;
        rolling = ((rolling << 2) + b - 1) & mArraySize;
        assert rolling <= mArraySize;
      }

      if (fill >= mWordSize) {
        //        System.err.println(templatePos + " store -> " + rolling);
        if ((templatePos - zeroBasedStart >= 0) && (templatePos - zeroBasedStart < mWords.length)) {
          mWords[templatePos - zeroBasedStart] = rolling + 1;
          mOccurs[rolling]++;
        }
      }
      ++templatePos;
    }

    //    System.err.println("end of first; " + templatePos + "fill = " + fill);
    // now step through the read/template
    int templateLast = -Math.min(offset, read.length);
    int readPos = 0;
    int readrolling = 0;
    int readfill = 0;
    int lb = 0;
    while (readPos < rlen) {
      //      System.err.println("second off: " + templatePos + " fill=" + fill);
      byte b = 0;
      if (templatePos >= 0 && templatePos < template.length) {
        b = template[templatePos]; // if it's off the end of the array, treat as N
        if (b < 1 || b > 4) { // if it's an N or something weird reset
          //System.err.println("B" + templatePos + ":" + b + ":" + template.length);
          return 0;
        }
      }
      if (b < 1 || b > 4) { // if it's an N or something weird reset
        fill = 0;
      } else { // otherwise count the bit and stream the hash
        ++fill;
        rolling = ((rolling << 2) + b - 1) & mArraySize;
      }

      if (fill >= mWordSize) {
        //        System.err.println(templatePos + " store -> " + rolling);
        if (templatePos - zeroBasedStart < mWords.length) {
          mWords[templatePos - zeroBasedStart] = rolling + 1;
        }
        mOccurs[rolling]++;
      }

      final byte readb = read[readPos + rStartPos];
      if (readb < 1 || readb > 4) { // if it's an N or something weird reset
        //System.err.println("C");
        return 0;
        //        readfill = 0;
        //        readrolling = 0;
      } else { // otherwise count the bit and stream the hash
        ++readfill;
        readrolling = ((readrolling << 2) + readb - 1) & mArraySize;
      }

      if (readfill >= mWordSize) {
        //        System.err.println("looking for " + readPos + " " + readrolling);
        if (mOccurs[readrolling] == 0) {
          //final int rp = readPos & (mWordSize - 1);
          final int rp = readPos % mWordSize;
          //          System.err.println("increasing ++ " + templatePos +" " + readPos + " % " + rp);
          mCount[rp]++;
          if (mCount[rp] > lb) {
            lb = mCount[rp];
            if (lb * mSubstitutionPenalty > maxScore) {
              break;
            }
          }
        }
      }
      // remove the last count (to the far left)
      if ((fill >= mWordSize) && (templateLast >= 0)) {
        if (mWords[templateLast] > 0) {
          mOccurs[mWords[templateLast] - 1]--;

          if (mOccurs[mWords[templateLast] - 1] < 0) {
            throw new RuntimeException(readPos + " " + templateLast + " " + zeroBasedStart + " " + (mWords[templateLast] - 1) + " "
                + mOccurs[mWords[templateLast] - 1]);
            //System.err.println("BUG");
          }
        }
      }
      if (templateLast >= 0) {
        mWords[templateLast] = -1; // set the word to -1 in case we try and use it
      }
      if (templatePos < template.length) {
        ++templatePos;
      }
      ++templateLast;
      ++readPos;
    }
    //int i = 0;
    for (final int element : mCount) {
      //System.err.println(i++ + " " + element);
      if (element > lb) {
        lb = element;
      }
    }
    return lb * mSubstitutionPenalty;
  }

  /**
   * Lower bound estimator, using a 4-nt lock-step creeping count.
   *
   * @param read the read
   * @param rlen read length
   * @param template the template
   * @param zeroBasedStart start position
   * @param maxScore maximum edit distance
   * @param maxShift the maximum allowed shift in start and end positions
   * @return alignment score
   */
  protected int calcLB(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift) {
    return calcLB(read, 0, rlen, template, zeroBasedStart, maxScore, maxShift);
  }
}
