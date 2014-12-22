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
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Path;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass("com.rtg.assembler.graph.implementation.GraphImplementationTest")
class PathLocal extends IntegralAbstract implements Path {

  private final GraphImplementation mGraphImplementation;
  private final long mPath;

  /**
   * @param path signed path identifier.
   * @param graphImplementation the graph that is parent of this (was an inner class before being moved out)
   */
  PathLocal(GraphImplementation graphImplementation, long path) {
    mGraphImplementation = graphImplementation;
    mPath = path;
  }

  @Override
  public boolean integrity() {
    mGraphImplementation.absPath(mPath);
    return true;
  }

  @Override
  public int length() {
    return mGraphImplementation.pathLength(mPath);
  }

  @Override
  public long contig(int index) {
    return mGraphImplementation.pathContig(mPath, index);
  }

  @Override
  public int index(long contigId) {
    mGraphImplementation.absContig(contigId);
    for (int i = 0; i < length(); i++) {
      if (contig(i) == contigId) {
        return i;
      }
    }
    return -1;
  }

  @Override
  public String toString() {
    final List<Long> contigs = new ArrayList<>();
    for (int i = 0; i < length(); i++) {
      contigs.add(mGraphImplementation.pathContig(mPath, i));
    }
    return "PathLocal{" + contigs + '}';
  }
}
