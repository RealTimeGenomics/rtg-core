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
package com.rtg.reader;

import java.util.Arrays;

/**
 * Given a sequence id return the <code>file</code> it belongs to
 */
public abstract class PointerFileLookup {

  private final int[] mStartSeqs;

  /**
   * Constructor
   * @param pointers pointers from {@link DataInMemory}
   */
  private PointerFileLookup(int[][] pointers) {
    mStartSeqs = new int[pointers.length];
    int startSeq = 0;
    for (int i = 0; i < mStartSeqs.length; i++) {
      mStartSeqs[i] = startSeq;
      startSeq += pointers[i].length - 1;
    }
  }

  /**
   * @param sequenceId sequence id within {@link SequencesReader}
   * @return logical file number that sequence is contained within
   */
  public abstract int lookup(int sequenceId);

  /**
   * @param fileNo file number to enquire about
   * @return the first sequence id that starts within this file (or undefined for files within which no sequences start)
   */
  public int startSeq(int fileNo) {
    return mStartSeqs[fileNo];
  }

  /**
   * @param pointers pointers from {@link DataInMemory}
   * @return pointer lookup class
   */
  public static PointerFileLookup generateLookup(int[][] pointers) {
    if (pointers.length <= 1) {
      return new SimplePointerFileLookup(pointers.length == 0 ? -1 : Integer.MAX_VALUE, pointers);
    } else if (pointers.length == 2) {
      return new SimplePointerFileLookup(pointers[0].length - 1, pointers);
    } else {
      return new BinarySearchPointerFileLookup(pointers);
    }
  }

  private static class SimplePointerFileLookup extends PointerFileLookup {
    private final int mTilt;
    SimplePointerFileLookup(int tilt, int[][] pointers) {
      super(pointers);
      mTilt = tilt;
    }

    @Override
    public int lookup(int sequenceId) {
      return sequenceId < mTilt ? 0 : 1;
    }
  }

  private static final class BinarySearchPointerFileLookup extends PointerFileLookup {
    private final int[] mLookup;

    private BinarySearchPointerFileLookup(int[][] pointers) {
      super(pointers);
      mLookup = new int[pointers.length];
      int tot = 0;
      for (int i = 0; i < pointers.length; i++) {
        mLookup[i] = tot;
        tot += pointers[i].length - 1;
      }
    }

    @Override
    public int lookup(int sequenceId) {
      int i = Arrays.binarySearch(mLookup, sequenceId);
      if (i >= 0) {
        while (i < mLookup.length && mLookup[i] == sequenceId) {
          i++;
        }
        return i - 1;
      } else {
        //i = -insertionPoint - 1
        //insertionPoint = -i - 1
        //insertionPoint = -(i + 1) [to avoid overflow]
        //we want insertion point - 1
        //ret = -(i + 1) - 1
        //    = -(i + 2)
        return -(i + 2);
      }
    }
  }
}
