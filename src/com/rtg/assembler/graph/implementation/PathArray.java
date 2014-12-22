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
    for (int i = 0; i < mContigs.length; i++) {
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
    for (int i = 0; i < mContigs.length; i++) {
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
    final List<Long> result = new ArrayList<>();
    for (int i = 0; i < p.length(); i++) {
      result.add(p.contig(i));
    }
    return result;
  }
}
