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
    for (int i = 0; i < length(); ++i) {
      if (contig(i) == contigId) {
        return i;
      }
    }
    return -1;
  }

  @Override
  public String toString() {
    final List<Long> contigs = new ArrayList<>(length());
    for (int i = 0; i < length(); ++i) {
      contigs.add(mGraphImplementation.pathContig(mPath, i));
    }
    return "PathLocal{" + contigs + '}';
  }
}
