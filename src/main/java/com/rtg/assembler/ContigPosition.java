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
