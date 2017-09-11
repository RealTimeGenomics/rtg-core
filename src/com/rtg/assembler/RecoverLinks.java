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

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.ContigByte;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;

/**
 */
public final class RecoverLinks {
  /** Attribute marking contigs created by Recover Links */
  public static final String RECOVERED_ATTRIBUTE = "recoveredLink";
  /** Description of recovered attribute */
  public static final String RECOVERED_DESC = "this contig has been created by paired end link recovery";

  //NOTE: the below is to be used in GraphMapTask if ever instituted
//  static final boolean RECOVER_LINKS = true; //Boolean.parseBoolean(System.getProperty("com.rtg.assembler.RecoverLinks", "true"));

  private RecoverLinks() { }
  static void recover(List<ConstraintCache> constraints, MutableGraph graph) {
    if (!graph.contigAttributes().containsKey(RECOVERED_ATTRIBUTE)) {
      graph.addContigAttribute(RECOVERED_ATTRIBUTE, RECOVERED_DESC);
    }
    if (!graph.pathAttributes().containsKey(RECOVERED_ATTRIBUTE)) {
      graph.addPathAttribute(RECOVERED_ATTRIBUTE, RECOVERED_DESC);
    }
    final ConstraintCache combined = ConstraintCache.combineCaches(constraints);
    final Set<Long> resolved = new HashSet<>();
    for (Map.Entry<Long, List<ConstraintCollector>> entry : combined.mIndex.entrySet()) {
      final long contigId = entry.getKey();
      // Skip contigs we've already examined from the other side
      if (resolved.contains(contigId)) {
        continue;
      }
      final List<ConstraintCollector> list = entry.getValue();
      if (list.size() != 1) {
        continue;
      }
      final ConstraintCollector constraint = list.get(0);
      if (combined.mIndex.get(constraint.mContigB).size() != 1) {
        continue;
      }
//      System.err.println("adding : " + constraint.mContigA + " -> " + constraint.mContigB);
      final int insert = constraint.averageDistance();
      // TODO  attempt aligning overlapping ends
      // We'll represent the link by creating a new contig that contains an overlap with each side and at least 2 Ns
      final int nPadding = Math.max(2, insert);
      final int overlap = graph.contigOverlap();
      final byte[] contigBytes = new byte[overlap * 2 + nPadding];
      final Contig start = graph.contig(constraint.mContigA);
      final Contig end = graph.contig(-constraint.mContigB);
      final int startOverlap = start.length() - overlap;
      for (int i = 0; i < overlap; ++i) {
        contigBytes[i] = start.nt(startOverlap + i);
      }
      final int endOverlap = contigBytes.length - overlap;
      for (int i = 0; i < overlap; ++i) {
        contigBytes[endOverlap + i] = end.nt(i);
      }
      final long newContig = graph.addContig(new ContigByte(contigBytes));
      final long newPath = graph.addPath(new PathArray(constraint.mContigA, newContig, -constraint.mContigB));
      final String readCount = String.valueOf(constraint.mConstraint.size());
      if (graph.pathAttributes().containsKey(GraphKmerAttribute.READ_COUNT)) {
        graph.setPathAttribute(newPath, GraphKmerAttribute.READ_COUNT, readCount);
      }
      if (graph.contigAttributes().containsKey(GraphKmerAttribute.READ_COUNT)) {
        graph.setContigAttribute(newContig, GraphKmerAttribute.READ_COUNT, readCount);
      }
      graph.setContigAttribute(newContig, RECOVERED_ATTRIBUTE, "true");
      resolved.add(constraint.mContigB);
      resolved.add(constraint.mContigA);
    }

  }
}
