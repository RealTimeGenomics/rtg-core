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

package com.rtg.assembler;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import com.rtg.assembler.graph.Graph;
import com.rtg.index.Finder;
import com.rtg.index.IndexCompressed;
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
    for (long i = 1; i <= g.numberContigs(); i++) {
      if (!g.contigDeleted(i)) {
        totalHashes += g.contigLength(i) / stepSize;
        totalLength += g.contigLength(i);
      }
    }
    assert totalLength > 0;
    final int valueBits = (int) Math.ceil(Math.log(totalLength) / Math.log(2));
    final int hashBits = wordSize * 2;
    final CreateParams create =  new CreateParams(totalHashes, hashBits, hashBits, valueBits, true, true, false, false);
    mIndex = new IndexCompressed(create, Integer.MAX_VALUE, false, Integer.MAX_VALUE, 0, 1);
    mWordSize = wordSize;
    mBuildFunction = new ExactHashFunction(wordSize, 2);
    for (int pass = 0; pass < 2; pass++) {
      int soFar = 0;
      for (long i = 1; i <= g.numberContigs(); i++) {
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
    for (int i = 0; i < graph.contigLength(contig); i++) {
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
    public void found(long id) {
      mHits.add(id);
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
   * @return A list of lists where each sub-list contains the contigs hit by the hash ending at that read position
   * @throws IOException
   */
  List<List<ContigPosition>> hits(byte[] read, Graph graph, ExactHashFunction searchFunction) throws IOException {
    final List<List<ContigPosition>> positions = new ArrayList<>();
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
