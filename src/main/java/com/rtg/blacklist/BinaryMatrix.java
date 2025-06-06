/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.blacklist;

import java.util.Arrays;
import java.util.Random;

/**
 * Represents a random square matrix of binary values (0 or 1) up to 64 long.
 */
public class BinaryMatrix {
  private static final int INITIAL_SEED = 13;
  private static final int MAX_ITERATIONS = 500;
  private final long[] mMatrix;
  private long[] mMatrixT = null; // For fast multiplication
  private final int mK;

  /**
   * @param k desired number of rows and columns
   */
  public BinaryMatrix(int k) {
    this(k, new Random(INITIAL_SEED));
  }

  BinaryMatrix(int k, Random r) {
    mK = k;
    mMatrix = new long[k];
    mMatrixT = new long[k];
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < k; ++j) {
        if (r.nextBoolean()) {
          mMatrix[i] |= 1L << j;
          mMatrixT[j] |= 1L << i;
        }
      }
    }
  }

  private BinaryMatrix(long[] matrix) {
    mMatrix = matrix;
    mMatrixT = null;
    mK = matrix.length;
  }

  private void transpose() {
    if (mMatrixT == null) {
      mMatrixT = new long[mK];
      for (int i = 0; i < mK; ++i) {
        long v = 0;
        for (int j = mK - 1; j >= 0; --j) {
          v <<= 1;
          v |= get(j, i);
        }
        mMatrixT[i] = v;
      }
    }
  }

  /**
   * Creates a random invertible matrix.
   * @param k desired number of rows and columns
   * @return the matrix
   */
  public static BinaryMatrix createReversibleMatrix(int k) {
    int seed = INITIAL_SEED;
    while (seed - INITIAL_SEED < MAX_ITERATIONS) {
      final BinaryMatrix b = new BinaryMatrix(k, new Random(seed));
      if (b.invert() != null) {
        return b;
      }
      ++seed;
    }
    throw new RuntimeException("Unable to create invertible matrix for size: " + k);
  }

  /**
   * Multiply given vector by the matrix
   * @param val the vector to multiply (assumed to be size k)
   * @return the resulting vector
   */
  public long times(final long val) {
    // Direct multiplication but with bit twiddling
    //    long res = 0;
    //    for (int i = 0; i < mK; ++i) {
    //      final long tmp = mMatrix[i] & val;
    //      res |= (Long.bitCount(tmp) & 1L) << i;
    //    }
    //    return res;

    // Multiplying via the transpose to avoid the bitCount
    long res = 0;
    long u = val;
    for (final long r : mMatrixT) {
      // if ((u & 1) != 0) {
      //   res ^= r;
      // }
      res ^= r & -(u & 1); // Avoid conditional logic
      u >>>= 1;
    }
    return res;
  }

  private long get(int j, int k) {
    return (mMatrix[j] >>> k) & 1;
  }

  private void swapRows(int a, int b) {
    final long tmp = mMatrix[a];
    mMatrix[a] = mMatrix[b];
    mMatrix[b] = tmp;
  }


  private static BinaryMatrix identity(int k) {
    final long[] m = new long[k];
    long j = 1;
    for (int i = 0; i < k; ++i, j <<= 1) {
      m[i] = j;
    }
    return new BinaryMatrix(m);
  }

  /**
   * Create the inverse for this matrix
   * @return the inverse matrix if it exists, null otherwise
   */
  public BinaryMatrix invert() {
    final BinaryMatrix a = new BinaryMatrix(Arrays.copyOf(mMatrix, mMatrix.length));
    final BinaryMatrix b = identity(mK);
    for (int k = 0; k < mK; ++k) {
      // Find first non-zero entry in kth column
      int j;
      for (j = k; j < mK; ++j) {
        if (0 != a.get(j, k)) {
          break;
        }
      }
      if (j == mK) {
        // Entire columns was zero, implies no solution
        return null;
      }
      if (j != k) {
        a.swapRows(j, k);
        b.swapRows(j, k);
      }
      // Add multiples of kth row to lower rows to zero entries
      for (int i = k + 1; i < mK; ++i) {
        final long u = a.get(i, k);
        if (u != 0) {
          a.mMatrix[i] ^= a.mMatrix[k];
          b.mMatrix[i] ^= b.mMatrix[k];
        }
      }
    }
    // Now do the Jordan step
    for (int k = mK - 1; k > 0; --k) {
      for (int j = k - 1; j >= 0; --j) {
        final long u = a.get(j, k);
        if (u != 0) {
          a.mMatrix[j] ^= a.mMatrix[k];
          b.mMatrix[j] ^= b.mMatrix[k];
          // Add multiples of the kth row to rows above to introduce zeros
        }
      }
    }
    b.transpose();
    return b;
  }

  /**
   * @return string representation of this matrix
   */
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < mK; ++i) {
      final String binString = Long.toBinaryString(mMatrix[i]);
      final String paddedString = binString.length() < mK ? String.format("%0" + (mK - binString.length()) + "d%s", 0, binString) : binString;
      sb.append(paddedString).append("\n");
    }
    return sb.toString();
  }
}
