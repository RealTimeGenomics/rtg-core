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
package com.rtg.metagenomics.metasnp;

import java.util.Arrays;
import java.util.List;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
*/
class ComplicatedBeta implements ProbAlpha {
  final int[] mCounts;
  final int[] mStride;
  final long mTotal;
  final int mStrains;
  ComplicatedBeta(List<Integer> refBytes, List<int[]> assignments, long length) {
    mStrains = assignments.get(0).length;
    mStride = new int[mStrains];
    int possibilities = 1;
    for (int i = 2; i <= mStrains + 1; ++i) {
      possibilities *= Math.min(i, AlphaSelector.MAX_VALUE + 1);
    }
    int stride = 1;
    // TODO potentially invert this stride array so that memory access pattern improves
    for (int i = 0; i < mStrains; ++i) {
      stride *= Math.min(i + 1, AlphaSelector.MAX_VALUE + 1);
      mStride[i] = stride;
    }
    mCounts = new int[possibilities];
    //TODO This initialisation could be improved the laplace correction includes some entries that are inaccessible
    mTotal = length + init(mCounts, mStrains, mStride);

    for (int position = 0; position < assignments.size(); ++position) {
      final int[] current = assignments.get(position);
      final int ref = refBytes.get(position);
      final int index = findIndex(ref, current);
      mCounts[index]++;
    }
  }
  public static int init(int[] counts, int strains, int[] strides) {
    int totalCells = 0;
    final int[] assignment = new int[strains];
    int currentIndex = 0;
    int last = -1;
    while (currentIndex > 0 || last <= max(assignment, currentIndex)) {
      assignment[currentIndex++] = last + 1;
      while (currentIndex < assignment.length) {
        assignment[currentIndex++] = 0;
      }
      int countIndex = 0;
      ++totalCells;
      for (int i = 0; i < strides.length; ++i) {
        countIndex += assignment[i] * strides[i];
      }
      counts[countIndex]++;
      last = assignment[--currentIndex];
      while (last == max(assignment, currentIndex) + 1 && currentIndex > 0) {
        last = assignment[--currentIndex];
      }
    }
    return totalCells;
  }

  private int findIndex(int ref, int[] assignment) {
    final int[] currentId = new int[mStrains];
    int index = 0;
    for (int i = 0; i < assignment.length; ++i) {
      int id = 0;
      if (assignment[i] == ref) {
        id = -1;
      } else {
        for (int j = 0; j < i; ++j) {
          if (assignment[i] == assignment[j]) {
            id = currentId[j] - 1;
            break;
          }
          id = Math.max(id, currentId[j]);
        }
      }
      currentId[i] = id + 1;
      index += currentId[i] * mStride[i];
    }
    return index;
  }
  int getCount(int ref, int[] assignment) {

    return mCounts[findIndex(ref, assignment)];
  }

  /**
   * @return the beta for the an assignment
   */
  @Override
  public double pAlpha(int ref, int[] current) {
    return (double) mCounts[findIndex(ref, current)] / mTotal;
  }

  public static int max(int[] arr, int limit) {
    int max = 0;
    for (int i = 0; i < limit; ++i) {
      final int current = arr[i];
      max = Math.max(max, current);
    }
    return Math.min(max, 2);
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    final int[] assignment = new int[mStrains];
    int currentIndex = 0;
    int last = -1;
    while (currentIndex > 0 || last <= max(assignment, currentIndex)) {
      assignment[currentIndex++] = last + 1;
      while (currentIndex < assignment.length) {
        assignment[currentIndex++] = 0;
      }
      int countIndex = 0;
      for (int i = 0; i < mStride.length; ++i) {
        countIndex += assignment[i] * mStride[i];
      }
      sb.append(Arrays.toString(assignment));
      sb.append(" = ");
      sb.append(mCounts[countIndex]);
      sb.append("\t");
      sb.append(Utils.realFormat((double) mCounts[countIndex] / mTotal, 6));
      sb.append(StringUtils.LS);
      last = assignment[--currentIndex];
      while (last == max(assignment, currentIndex) + 1 && currentIndex > 0) {
        last = assignment[--currentIndex];
      }
    }
    return sb.toString();
  }
}
