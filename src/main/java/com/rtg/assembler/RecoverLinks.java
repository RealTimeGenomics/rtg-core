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
