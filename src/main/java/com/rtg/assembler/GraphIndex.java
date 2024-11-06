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

package com.rtg.assembler;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import com.rtg.assembler.graph.Graph;
import com.rtg.index.Finder;
import com.rtg.index.IndexCompressed;
import com.rtg.index.UnfilteredFilterMethod;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.params.CreateParams;

/**
 */
public class GraphIndex {
  final IndexCompressed mIndex;
  final ExactHashFunction mBuildFunction;
  final int mWordSize;
  final TreeMap<Long, Long> mPositionDecoder;

  GraphIndex(Graph g, int stepSize, int wordSize) {
    mPositionDecoder = ContigPosition.buildDecoder(g);
    long totalHashes = 0;
    long totalLength = 0;
    for (long i = 1; i <= g.numberContigs(); ++i) {
      if (!g.contigDeleted(i)) {
        totalHashes += g.contigLength(i) / stepSize;
        totalLength += g.contigLength(i);
      }
    }
    assert totalLength > 0;
    final int valueBits = (int) Math.ceil(Math.log(totalLength) / Math.log(2));
    final int hashBits = wordSize * 2;
    final CreateParams create =  new CreateParams(totalHashes, hashBits, hashBits, valueBits, true, true, false, false);
    mIndex = new IndexCompressed(create, new UnfilteredFilterMethod(), 1);
    mWordSize = wordSize;
    mBuildFunction = new ExactHashFunction(wordSize, 2);
    for (int pass = 0; pass < 2; ++pass) {
      int soFar = 0;
      for (long i = 1; i <= g.numberContigs(); ++i) {
        if (!g.contigDeleted(i)) {
          addHashes(g, i, stepSize, soFar);
          soFar +=  g.contigLength(i);
        }

      }
      mIndex.freeze();
    }
  }

  ExactHashFunction getSearchFunction() {
    return new ExactHashFunction(mWordSize, 2, true);
  }

  /**
   * @return the size of the window use in searching
   */
  public int windowSize() {
    return mWordSize;
  }

  private void addHashes(Graph graph, long contig, int stepSize, int soFar) {
    mBuildFunction.reset();
    for (int i = 0; i < graph.contigLength(contig); ++i) {
      final byte code = (byte) (graph.nt(contig, i) - 1);
      if (code < 0) {
        mBuildFunction.reset();
      } else {
        mBuildFunction.hashStep(code);
      }
      if (mBuildFunction.isValid() && i % stepSize == (stepSize  - 1)) {
        mIndex.add(mBuildFunction.hash(), soFar + i);
      }
    }
  }
  static class HitList extends Finder {
    List<Long> mHits = new ArrayList<>();

    @Override
    public boolean found(long id) {
      mHits.add(id);
      return true;
    }
  }

  List<Long> hitContigs(long hash) throws IOException {
    final HitList finder = new HitList();
    mIndex.search(hash, finder);
    return finder.mHits;
  }

  List<ContigPosition> hashHits(ExactHashFunction hash, Graph graph) throws IOException {
    if (!hash.isValid()) {
      return Collections.emptyList();
    }
    final GraphIndex.HitList hitListForward = new GraphIndex.HitList();
    final GraphIndex.HitList hitListReverse = new GraphIndex.HitList();
    mIndex.search(hash.hash(), hitListForward);
    mIndex.search(hash.hashReverse(), hitListReverse);
    final List<ContigPosition> hitList = new ArrayList<>();
    for (long hit : hitListForward.mHits) {
      hitList.add(new ContigPosition(hit, graph, mPositionDecoder));
    }
    for (long hit : hitListReverse.mHits) {
      final ContigPosition forward = new ContigPosition(hit, graph, mPositionDecoder);
      // reverse the contig
      int position = graph.contigLength(forward.mContigId) - forward.mPosition - 1;
      // move the hit position to the other end of the hash
      position = position + hash.getWindowSize() - 1;
      hitList.add(new ContigPosition(-forward.mContigId, position, graph));
    }
    return hitList;

  }

  /**
   * Collect all the hits that occur for each window in this read
   * @param read byte representation of the read
   * @param graph the graph to search
   * @param searchFunction the hashing function to use for searching
   * @return A list of lists where each sub-list contains the contigs hit by the hash ending at that read position
   * @throws IOException if an I/O error occurs
   */
  List<List<ContigPosition>> hits(byte[] read, Graph graph, ExactHashFunction searchFunction) throws IOException {
    final List<List<ContigPosition>> positions = new ArrayList<>(read.length);
    searchFunction.reset();
    for (byte base : read) {
      final byte code = (byte) (base - 1);
      if (code < 0) {
        searchFunction.reset();
      } else {
        searchFunction.hashStep(code);
      }
      positions.add(hashHits(searchFunction, graph));
    }
    assert positions.size() == read.length;
    return positions;
  }
}
