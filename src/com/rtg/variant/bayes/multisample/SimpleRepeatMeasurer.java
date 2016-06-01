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
package com.rtg.variant.bayes.multisample;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;

/**
 * Measures the total length of simple repeat regions between two positions (may consider nucleotides outside this range to determine repeat length).
 * Simple repeats are anything where a string up to 3 length occurs over and over.
 *
 */
public class SimpleRepeatMeasurer implements RepeatMeasurer {

  static final int MAX_MER_LENGTH = GlobalFlags.getIntegerValue(CoreGlobalFlags.COMPLEX_REGION_SIMPLE_REPEAT_LIMIT);

  protected final byte[] mReferenceNts;

  protected final int mMaxMerLength;

  /**
   * Constructor
   * @param referenceNts reference sequence bytes
   */
  public SimpleRepeatMeasurer(byte[] referenceNts) {
    this(referenceNts, MAX_MER_LENGTH);
  }

  /**
   * Constructor
   * @param referenceNts reference sequence bytes
   * @param maxMerLength length of longest n-mer simple repeat type
   */
  public SimpleRepeatMeasurer(byte[] referenceNts, int maxMerLength) {
    mReferenceNts = referenceNts;
    mMaxMerLength = maxMerLength;
  }

  /**
   * Measures the total length of simple repeat regions between two positions (may consider nucleotides outside this range to determine repeat length).
   * Simple repeats are anything where a string up to 3 length occurs over and over.
   *
   * <code>accaccac</code> will count as a simple repeat of length 8; The third iteration is only partial but still counts.
   * <code>accac</code> <code>5nt</code>
   * <code>acca</code> is <code>2nt</code> of repeat, unless followed by <code>c</code>, in which case it's <code>4nt</code> of repeat
   * <code>aaa</code> is <code>3nt</code> of repeat
   * <code>aca</code> is not any repeat, unless followed by <code>c</code>, in which case it's <code>3nt</code> of repeat
   * <code>acc</code> is <code>2nt</code> of repeat (cc)
   * <code>tt</code> is <code>2nt</code> of repeat
   * <code>at</code> is not any repeat.
   * <code>aatcagtt</code> is <code>2+2nt</code> of repeat (aa and <code>tt</code>).
   *
   * @param positionA template start position
   * @param positionB template end position (exclusive)
   * @return the total length of repeats
   */
  @Override
  public int measureRepeats(int positionA, int positionB) {
    final int clippedPositionA = Math.max(positionA, 0);
    int repeatTotal = 0;
    //StringBuilder sb = new StringBuilder();
    for (int repeatStart = clippedPositionA; repeatStart < positionB; repeatStart++) {
      for (int simpleLength = 1; simpleLength <= mMaxMerLength; simpleLength++) {
        // work out the potential repeat unit
        if (mReferenceNts.length > repeatStart + simpleLength) {
          int repeatPosition = repeatStart + simpleLength;
          if (simpleLength < 3) { //1 & 2 long repeats
            while (repeatPosition < positionB && isRepeat(repeatStart, simpleLength, repeatPosition)) {
              repeatPosition += simpleLength;
            }
            if (repeatPosition != repeatStart + simpleLength) {
              final int repeatEnd = Math.min(repeatPosition + mismatchPosition(repeatStart, simpleLength, repeatPosition), positionB);
              //sb.append("r" + simpleLength + ':' + repeatStart + '-' + repeatEnd + " ");
              repeatTotal += repeatEnd - repeatStart;
              repeatStart = repeatPosition + mismatchPosition(repeatStart, simpleLength, repeatPosition);
              break;
            }
          } else {  //special case for detecting 3 long repeats
            int mmlength = 0;
            while (repeatPosition < positionB && (mmlength = mismatchPosition(repeatStart, simpleLength, repeatPosition)) == simpleLength) { // isRepeat(repeatStart, simpleLength, repeatPosition)) {
              repeatPosition += simpleLength;
              mmlength = 0;
            }
            if (repeatPosition != repeatStart + simpleLength || mmlength >= 2) {
              final int repeatEnd = Math.min(repeatPosition + mismatchPosition(repeatStart, simpleLength, repeatPosition), positionB);
              //sb.append("r" + simpleLength + ':' + repeatStart + '-' + repeatEnd + " ");
              repeatTotal += repeatEnd - repeatStart;
              repeatStart = repeatPosition + mismatchPosition(repeatStart, simpleLength, repeatPosition);
            }
          }
        }
      }
    }
    //System.err.println("Total repeats: " + repeatTotal + " via " + sb.toString());
    return repeatTotal;
  }

  /**
   * @param repeatStart start position of the first repeat
   * @param simpleLength length of the repeat
   * @param repeatPosition the start position of the repeat to check for mismatch
   * @return true if a repeat of <code>simpleLength</code> is found at the <code>repeatPosition</code>
   */
  protected boolean isRepeat(int repeatStart, int simpleLength, int repeatPosition) {
    return mismatchPosition(repeatStart, simpleLength, repeatPosition) == simpleLength;
  }

  /**
   * @param repeatStart start position of the first repeat
   * @param simpleLength length of the repeat
   * @param repeatPosition the start position of the repeat to check for mismatch
   * @return the position where the repeat is broken (or <code>simpleLength</code> if the end is reached)
   */
  protected int mismatchPosition(int repeatStart, int simpleLength, int repeatPosition) {
    int i;
    for (i = 0; i < simpleLength && repeatPosition + i < mReferenceNts.length; i++) {
      if (mReferenceNts[repeatStart + i] != mReferenceNts[repeatPosition + i]) {
        return i;
      }
    }
    return i;
  }

}
