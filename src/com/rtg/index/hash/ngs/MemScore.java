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
package com.rtg.index.hash.ngs;

import com.rtg.mode.DNA;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Score given potential alignment given an anchor position.
 *
 */
public class MemScore extends IntegralAbstract {

  /**
   * Count the number of mismatches in a putative alignment.  That is, the lower
   * the return value, the better the match.  Positions are counted from 1.
   *
   * @param template template sequence data
   * @param tposn start anchor position in template
   * @param read read sequence data
   * @param readposn start anchor position in read
   * @param reverseComplement true iff match is reverse complement
   * @return a score
   */
  public static int score(final byte[] template, final int tposn, final byte[] read, final int readposn, final boolean reverseComplement) {
    assert tposn > 0 && tposn <= template.length;
    assert readposn > 0 && readposn <= read.length : readposn;

    final int toffset = tposn - 1;
    final int roffset = readposn - 1;

    // todo: special case where roffset=0 && length==read.length, more generally
    // we should not need to recheck values already known to agree, i.e. within
    // the passed length
    // todo: handle simple indels

    int mismatchCount = 0;

    if (reverseComplement) {
      // match rightwards along read
      for (int k = roffset, j = toffset; k < read.length; k++, j--) {
        if (j < 0 || template[j] != DNA.complement(read[k])) {
          mismatchCount++;
        }
      }
      // match leftwards along read
      for (int k = roffset - 1, j = toffset + 1; k >= 0; k--, j++) {
        if (j >= template.length || template[j] != DNA.complement(read[k])) {
          mismatchCount++;
        }
      }
    } else {
      // match rightwards along read
      for (int k = roffset, j = toffset; k < read.length; k++, j++) {
        if (j >= template.length || template[j] != read[k]) {
          mismatchCount++;
        }
      }
      // match leftwards along read
      for (int k = roffset - 1, j = toffset - 1; k >= 0; k--, j--) {
        if (j < 0 || template[j] != read[k]) {
          mismatchCount++;
        }
      }
    }
    return mismatchCount;
  }


  /**
   * @return bitwise difference in the arguments.
   */
  static long simpleDiff(final long a0, final long a1, final long b0, final long b1) {
    //System.err.println("a0=" + com.rtg.util.Utils.toBitsSep(a0));
    //System.err.println("a1=" + com.rtg.util.Utils.toBitsSep(a1));
    //System.err.println("b0=" + com.rtg.util.Utils.toBitsSep(b0));
    //System.err.println("b1=" + com.rtg.util.Utils.toBitsSep(b1));
    return (a0 ^ b0) | (a1 ^ b1);
  }

  private final int mReadLength;

  private final long mMask;

  /**
   * @param readLength length of the reads.
   */
  public MemScore(final int readLength) {
    mReadLength = readLength;
    mMask = readLength == 64 ? -1L : (1L << mReadLength) - 1;
  }

  /**
   * Compute a fast and possibly inaccurate score for the mismatch. (0 is good).
   * Accurate for any number of substitutions - does not deal with indels.
   * @param r0 bit 0 of reads
   * @param r1 bit 1 of reads
   * @param template0 high order bits of the recent template (may need to be masked as it contains all the recent codes up to a length of 64).
   * @param template1 low order bits of the recent template (may need to be masked as it contains all the recent codes up to a length of 64).
   * @return the score.
   */
  public int fastScore(final long r0, long r1, final long template0, final long template1) {
    //System.err.println("fastScore readId=" + readId);
    //System.err.println("t0=" + Utils.toBitsSep(template0));
    //System.err.println("t1=" + Utils.toBitsSep(template1));
    //System.err.println("ma=" + Utils.toBitsSep(mMask));

    final long diff = simpleDiff(r0, r1, template0, template1) & mMask;
    return Long.bitCount(diff);
  }

  /**
   * Compute a fast and possibly inaccurate score for the mismatch. (0 is good).
   * Accurate for a single indel, more indels and substitutions will be an approximation.
   * @param r0 bit 0 of reads
   * @param r1 bit 1 of reads
   * @param template0 high order bits of the recent template (may need to be masked as it contains all the recent codes up to a length of 64).
   * @param template1 low order bits of the recent template (may need to be masked as it contains all the recent codes up to a length of 64).
   * @return the score.
   */
  public int indelScore(final long r0, long r1, final long template0, final long template1) {
    //System.err.println("r0=" + Utils.toBitsSep(r0));
    //System.err.println("r1=" + Utils.toBitsSep(r1));
    final long diff = simpleDiff(r0, r1, template0, template1) & mMask;
    final int exactScore = Long.bitCount(diff);
    if (exactScore <= 2) {
      return exactScore;
    }
    final long diffDel = diff & simpleDiff(r0 >>> 1, r1 >>> 1, template0, template1);
    final int scoreDel = Long.bitCount(diffDel);
    final long diffIns = diff & simpleDiff(r0, r1, template0 >>> 1, template1 >>> 1);
    final int scoreIns = Long.bitCount(diffIns);
    //    System.err.println("r0=" + com.rtg.util.Utils.toBits(r0, mReadLength));
//    System.err.println("r1=" + com.rtg.util.Utils.toBits(r1, mReadLength));
//    System.err.println("t0=" + com.rtg.util.Utils.toBits(template0, mReadLength));
//    System.err.println("t1=" + com.rtg.util.Utils.toBits(template1, mReadLength));
//    System.err.println("d =" + com.rtg.util.Utils.toBits(diff, mReadLength) + " ex=" + exactScore);
//    System.err.println("dd=" + com.rtg.util.Utils.toBits(diffDel, mReadLength) + " sd=" + scoreDel);
//    System.err.println("di=" + com.rtg.util.Utils.toBits(diffIns, mReadLength) + " si=" + scoreIns);
//    System.err.println("be=" + bestScore);
    return Math.min(exactScore, Math.max(1, Math.min(scoreDel, scoreIns)) + 1);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mReadLength >= 1 && mReadLength <= 64);
    for (int i = 0; i < mReadLength; i++) {
      Exam.assertTrue("i=" + i + " readLength=" + mReadLength, (mMask & (1L << i)) != 0);
    }
    for (int i = mReadLength; i < 64; i++) {
      Exam.assertTrue((mMask & (1L << i)) == 0);
    }
    return true;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Score length=").append(mReadLength);
  }
}

