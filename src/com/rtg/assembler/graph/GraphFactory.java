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

package com.rtg.assembler.graph;

import java.util.Map;

import com.rtg.assembler.graph.implementation.GraphImplementation;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;

/**
 */
public interface GraphFactory {

  /**
   * Factory that uses the default graph implementation.
   */
  GraphFactory DEFAULT = new GraphFactory() {
    @Override
    public MutableGraph makeGraph(int contigOverlap, Map<String, String> contigAttributes, Map<String, String> pathAttributes) {
      return new GraphImplementation(contigOverlap, contigAttributes, pathAttributes);
    }
  };

  /**
   * Factory that constructs a graph with efficient support for the <code>kMerFreq</code> attribute
   */
  GraphFactory KMER = new GraphFactory() {
    @Override
    public MutableGraph makeGraph(int contigOverlap, Map<String, String> contigAttributes, Map<String, String> pathAttributes) {
      return new GraphKmerAttribute(contigOverlap, contigAttributes, pathAttributes);
    }
  };

  /**
   * Construct a new MutableGraph given the contig and path attributes.
   * @param contigOverlap number of bases adjacent contigs overlap by
   * @param contigAttributes contig attributes map from attribute name to description
   * @param pathAttributes path attributes map from attribute name to description
   * @return a new graph.
   */
  MutableGraph makeGraph(int contigOverlap, Map<String, String> contigAttributes, Map<String, String> pathAttributes);

}
