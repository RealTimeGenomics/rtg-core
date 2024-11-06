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

import java.util.ArrayList;
import java.util.List;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.MutableGraph;

/**
 * generates a contig as the concatenation of a series of contigs which overlap by <code>kmerSize - 1</code>
 * This will not create any additional objects.
 */
public class ContigConcatenated implements Contig {

  private final List<Contig> mContigs;
  private final int mKmerSize;
  private List<Long> mStartLinks;
  private List<Long> mEndLinks;

  /**
   * Construct a new contig
   * @param walkedIds the existing contig ids to include
   * @param graph graph containing the existing contigs
   * @param kmerSize the kmer size used to construct the graph
   */
  public ContigConcatenated(List<Long> walkedIds, MutableGraph graph, int kmerSize) {
    mContigs = new ArrayList<>(walkedIds.size());
    for (long id : walkedIds) {
      mContigs.add(graph.contig(id));
    }
    mKmerSize = kmerSize;
  }

  @Override
  public int length() {
    int totalLength = mKmerSize - 1;
    for (Contig contig : mContigs) {
      totalLength += contig.length() - mKmerSize + 1;
    }
    return totalLength;
  }

  @Override
  public byte nt(int index) {
    int soFar = 0;
    int overlap = 0;
    for (Contig contig : mContigs) {
      final int contigLength = contig.length();
      if (soFar + contigLength - overlap > index) {
       return contig.nt(index - soFar + overlap);
      }
      soFar += contigLength - overlap;
      overlap = mKmerSize - 1;
    }
    throw new IllegalArgumentException("requested an nt > length");
  }

  /**
   * @return the contigs linking to the start of this one
   */
  public List<Long> getStartLinks() {
    return mStartLinks;
  }

  /**
   * @param startLinks the contigs linking to the start of this one
   */
  public void setStartLinks(List<Long> startLinks) {
    mStartLinks = startLinks;
  }

  /**
   * @return the contigs linking to the end of this one
   */
  public List<Long> getEndLinks() {
    return mEndLinks;
  }

  /**
   * @param endLinks the contigs linking to the end of this one
   */
  public void setEndLinks(List<Long> endLinks) {
    mEndLinks = endLinks;
  }
}
