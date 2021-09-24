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
package com.rtg.index.hash.ngs.general;

import java.util.ArrayList;
import java.util.Collection;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Does the tricky calculations necessary for finding a mask which fits the
 * requested constraints.
 * To deal with the really tricky cases the effective length of the read may be reduced.
 */
public final class Skeleton extends IntegralAbstract {

  private final int mReadLength;

  private final int mWindowLength;

  private final int mSubstitutions;

  private final int mIndels;

  private final int mIndelLength;

  private final int mWindowLengthActual;

  private final int mReadLengthActual;

  private final int mChunkLength;

  private final int mChunks;

  private final int mWindowChunks;

  private final boolean mValid;

  /**
   * @param readLength lengths of reads.
   * @param windowLength minimum length of the window used for indexing.
   * @param substitutions minimum number of substitution that will be found.
   * @param indels minimum number of indels that will be found.
   * @param indelLength indel length
   */
  public Skeleton(final int readLength, final int windowLength, final int substitutions, final int indels, final int indelLength) {
    mIndels = indels;
    mIndelLength = indelLength;
    mReadLength = readLength;
    mSubstitutions = substitutions;
    mWindowLength = windowLength;
    final boolean inOk =
      mReadLength > 0 && mReadLength <= 64
      && mWindowLength > 0
      && mWindowLength <= mReadLength
      && mWindowLength + mSubstitutions <= 32
      && mSubstitutions >= 0 && mSubstitutions <= (mReadLength - mWindowLength)
      && mIndels >= 0 && mIndels <= mSubstitutions
      && mIndelLength >= 1
      ;
      if (!inOk) {
        mWindowLengthActual = -1;
        mReadLengthActual = -1;
        mChunkLength = -1;
        mChunks = -1;
        mWindowChunks = -1;
        mValid = false;
        return;
      }
      int k = 1;
      while (true) {
        if (k > mWindowLength) {
          throw new RuntimeException();
        }

        final int t = (mWindowLength + k - 1) / k;
        assert t >= 1;
        final int c = k + mSubstitutions;
        final int w = k * t;
        final int r = c * t;
        if (w >= mWindowLength && w <= 32 && r <= mReadLength) {
          mWindowLengthActual = w;
          mReadLengthActual = r;
          mChunkLength = t;
          mChunks = c;
          mWindowChunks = k;
          mValid = true;
          break;
        }
        ++k;
      }

      integrity();
  }

  /**
   * The total number of windows.
   * @return the number of windows.
   */
  public int numberWindows() {
    return (int) Util.binomial(mChunks, mSubstitutions);
  }

  /**
   * Check if inputs used for construction lead to a valid
   * skeleton.
   * @return true iff the constraints have been satisfied.
   */
  public boolean valid() {
    return mValid;
  }

  /**
   * Compute the set of masks (determined by the skeleton).
   * @return the set of masks.
   */
  public Collection<SingleMask> masks() {
    final Collection<SingleMask> masks = new ArrayList<>();
    //System.err.println("chunks=" + mChunks + " windowChunks=" + mWindowChunks);
    final Combinator comb = new Combinator(mChunks, mWindowChunks) {
      @Override
      public void permutation() {
        // create new SingleMask and populate it
        final SingleMask sm = new SingleMask(mChunks, mWindowLengthActual, mIndels, mIndelLength);
        boolean active = false;
        int start = 0;
        int count = 0;
        for (int i = 0; true; ++i) {
          final boolean set = i != mChunks && getBit(i);
          if (set) {
            ++count;
          }
          if (set && !active) {
            active = true;
            start = i;
          }
          if (!set && active) {
            active = false;
            final Skel sk = new Skel(i * mChunkLength - 1, (i - start) * mChunkLength, count * mChunkLength - 1);
            sm.add(sk);
          }
          if (i == mChunks) {
            break;
          }
        }
        assert count == mWindowChunks;
        sm.freeze();
        masks.add(sm);
      }
    };
    comb.combine();
    assert numberWindows() == masks.size();
    return masks;
  }

  /**
   * Get the length of the reads.
   * @return the length of the reads.
   */
  public int readLength() {
    return mReadLength;
  }

  /**
   * Get the number of nt in a window.
   * @return the number of nt in a window.
   */
  public int windowLength() {
    return mWindowLengthActual;
  }

  /**
   * Get the number of bits needed to uniquely represent a window (the hash may lose information).
   * @return the number of bits needed to uniquely represent a window.
   */
  public int windowBits() {
    return 2 * mWindowLengthActual;
  }

  /**
   * Get the number of substitutions that are guaranteed to be found.
   * @return the number of substitutions.
   */
  public int substitutions() {
    return mSubstitutions;
  }

  /**
   * Get the number of indels that are guaranteed to be found.
   * @return the number of indels.
   */
  public int indels() {
    return mIndels;
  }

  /**
   * Get the maximum length of an indel.
   * @return the maximum length of an indel.
   */
  int indelLength() {
    return mIndelLength;
  }

  /**
   * Get the length of individual chunks.
   * @return length of individual chunks.
   */
  int chunkLength() {
    return mChunkLength;
  }

  /**
   * Display a human readable layout of the details of the masks.
   * @return the layout.
   */
  public String dumpMask() {
    final StringBuilder sb = new StringBuilder();
    sb.append("" + "Mask set for readlength = ").append(mReadLength).append(" wordsize = ").append(mWindowLengthActual).append(" substitutions = ").append(mSubstitutions).append(" indels = ").append(mIndels).append(" indelLength = ").append(mIndelLength).append(LS).append(LS);

    for (final SingleMask mask : masks()) {
      long v = 0;
      for (int i = 0; i < mask.size(); ++i) {
        final Skel sk = mask.subSkeleton(i);
        if (sk != null) {
          final long bits = ((1L << sk.length()) - 1) << (sk.position() - sk.length() + 1);
          v |= bits;
        }
      }
      final String s = "0000000000000000000000000000000000000000000000000000000000000000" + Long.toBinaryString(v);
      sb.append(s.substring(s.length() - mReadLength)).append(LS);
    }
    return sb.toString();
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Inputs:").append(mValid ? "" : "invalid").append(LS);
    sb.append("  Read length   ").append(mReadLength).append(LS);
    sb.append("  Window length ").append(mWindowLength).append(LS);
    sb.append("  Substitutions ").append(mSubstitutions).append(LS);
    sb.append("  Indels        ").append(mIndels).append(LS);
    sb.append("  Indel length  ").append(mIndelLength).append(LS);
    if (!mValid) {
      return;
    }
    sb.append("Computed:").append(LS);
    sb.append("  Chunk length  ").append(mChunkLength).append(LS);
    sb.append("  Window chunks ").append(mWindowChunks).append(LS);
    sb.append("  Window actual ").append(mWindowLengthActual).append(LS);
    sb.append("  Read actual   ").append(mReadLengthActual).append(LS);
    sb.append("  Number masks  ").append(numberWindows()).append(LS);
  }

  @Override
  public boolean integrity() {
    if (!mValid) {
      Exam.assertEquals(mWindowLengthActual, -1);
      Exam.assertEquals(mReadLengthActual, -1);
      Exam.assertEquals(mChunkLength, -1);
      Exam.assertEquals(mChunks, -1);
      Exam.assertEquals(mWindowChunks, -1);
      return true;
    }
    Exam.assertTrue(mReadLength > 0 && mReadLength <= 64);
    Exam.assertTrue(mWindowLength > 0 && mWindowLength <= 32);
    Exam.assertTrue(mWindowLength <= mReadLength);
    Exam.assertTrue(mSubstitutions >= 0 && mSubstitutions <= (mReadLength - mWindowLength));
    Exam.assertTrue(mWindowLength + mSubstitutions <= 32);
    Exam.assertTrue(mIndels >= 0 && mIndels <= mSubstitutions);
    Exam.assertTrue(mIndelLength >= 1);

    Exam.assertTrue(mWindowLengthActual + this.toString(), mWindowLengthActual > 0 && mWindowLengthActual <= 32);
    Exam.assertTrue(mWindowLengthActual <= mReadLength);
    Exam.assertTrue(mWindowLengthActual >= mWindowLength);
    Exam.assertTrue(mChunkLength > 0);
    Exam.assertTrue(mReadLengthActual  <= mReadLength);
    Exam.assertTrue(this.toString(), mReadLengthActual == mChunks * mChunkLength);

    Exam.assertEquals(mWindowLengthActual, mWindowChunks * mChunkLength);
    Exam.assertTrue(mWindowLengthActual <= 32);
    return true;
  }

  /**
   * Handy way of getting a quick look at the skeleton that will be built.
   * @param args command line arguments: read length, window size, substitutions, indels.
   */
  public static void main(final String[] args) {
    final int r = Integer.parseInt(args[0]);
    final int w = Integer.parseInt(args[1]);
    final int s = Integer.parseInt(args[2]);
    final int i = Integer.parseInt(args[3]);
    final int l = Integer.parseInt(args[4]);
    final Skeleton sk = new Skeleton(r, w, s, i, l);
    sk.integrity();
    System.out.println(sk);
    final Collection<SingleMask> masks = sk.masks();
    int j = 0;
    for (SingleMask sm : masks) {
      System.out.println("[" + j + "]");
      System.out.println(sm);
      System.out.println();
      ++j;
    }
    System.out.println(sk.dumpMask());
  }
}
