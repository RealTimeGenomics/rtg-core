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
