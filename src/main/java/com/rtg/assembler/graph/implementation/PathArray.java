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

package com.rtg.assembler.graph.implementation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.assembler.graph.Path;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Intended as the implementation used when putting into a graph.
 */
public class PathArray extends IntegralAbstract implements Path {

  private final long[] mContigs;

  /**
   * @param contigs the contigs in the path.
   */
  public PathArray(long... contigs) {
    mContigs = contigs;
  }

  /**
   * Construct from the contig ids in a list
   * @param contigs the contigs in the path
   */
  public PathArray(List<Long> contigs) {
    mContigs = new long[contigs.size()];
    for (int i = 0; i < mContigs.length; ++i) {
      mContigs[i] = contigs.get(i);
    }
  }

  @Override
  public int length() {
    return mContigs.length;
  }

  @Override
  public long contig(int index) {
    return mContigs[index];
  }

  @Override
  public int index(long contigId) {
    if (contigId == 0) {
      throw new IllegalArgumentException();
    }
    for (int i = 0; i < mContigs.length; ++i) {
      if (mContigs[i] == contigId) {
        return i;
      }
    }
    return -1;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (final long contig : mContigs) {
      Exam.assertTrue(contig != 0);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mContigs.length >= 2);
    return true;
  }

  @Override
  public String toString() {
    return "PathArray{" + "mContigs=" + Arrays.toString(mContigs) + '}';
  }

  /**
   * Construct a list of the contigs in this path
   * @param p the path to process
   * @return a list containing all the contig ids in p in order
   */
  public static List<Long> toList(Path p) {
    final List<Long> result = new ArrayList<>(p.length());
    for (int i = 0; i < p.length(); ++i) {
      result.add(p.contig(i));
    }
    return result;
  }
}
