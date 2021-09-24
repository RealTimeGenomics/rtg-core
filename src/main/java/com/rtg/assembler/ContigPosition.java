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

import java.util.Map;
import java.util.TreeMap;

import com.rtg.assembler.graph.Graph;

/**
*/
class ContigPosition {

  static TreeMap<Long, Long> buildDecoder(Graph graph) {
    final TreeMap<Long, Long> positionDecoder = new TreeMap<>();
    long total = 0;
    for (long i = 1; i <= graph.numberContigs(); ++i) {
      if (!graph.contigDeleted(i)) {
        positionDecoder.put(total, i);
        total += graph.contigLength(i);
      }
    }
    return positionDecoder;
  }

  final long mContigId;
  final int mPosition;
  final Graph mGraph;

  ContigPosition(long contigId, int position, Graph graph) {
    mContigId = contigId;
    mPosition = position;
    mGraph = graph;
  }

  ContigPosition(long encoded,  Graph graph, TreeMap<Long, Long> positionDecoder) {
    final Map.Entry<Long, Long> e = positionDecoder.floorEntry(encoded);
    mContigId = e.getValue();
    mPosition = (int) (encoded - e.getKey());
    mGraph = graph;
  }

  long encode() {
    long total = 0;
    for (long i = 1; i < mContigId; ++i) {
      if (!mGraph.contigDeleted(i)) {
        total += mGraph.contigLength(i);
      }
    }
    return total + mPosition;
  }

  @Override
  public String toString() {
    return "ContigPosition{" + "mContigId=" + mContigId + ", mPosition=" + mPosition + '}';
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }

    final ContigPosition that = (ContigPosition) o;

    if (mContigId != that.mContigId) {
      return false;
    }
    if (mPosition != that.mPosition) {
      return false;
    }

    return true;
  }

  @Override
  public int hashCode() {
    int result = (int) (mContigId ^ (mContigId >>> 32));
    result = 31 * result + mPosition;
    return result;
  }
}
