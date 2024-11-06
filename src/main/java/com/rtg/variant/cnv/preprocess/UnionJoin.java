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
package com.rtg.variant.cnv.preprocess;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import com.rtg.util.MultiMap;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusComparator;

/**
 * Joins two datasets using the union of their regions inserting empty cells where necessary.
 * The joined datasets must be sorted.
 */
public class UnionJoin extends SimpleJoin {

  UnionJoin(File in) throws IOException {
    super(in);
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   */
  public UnionJoin(RegionDataset dataset, String prefix) {
    super(dataset, prefix);
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    final RegionDataset src = mDataset.copy();
    final List<SequenceNameLocus> union = unionRegions(src, dataset);
    expandRegions(union, src);
    expandRegions(union, dataset);
    simpleJoin(src, dataset);
  }

  // Get the list of all regions in both datasets.
  static List<SequenceNameLocus> unionRegions(RegionDataset d1, RegionDataset d2) {
    // This is quite yucky - since there isn't a reliable global ordering on sequence names,
    // we get the initial ordering from the first dataset and then insert missing entries
    // from the second dataset where it looks right.
    final List<String> chrs = new ArrayList<>();
    final MultiMap<String, SequenceNameLocus> unionByChr = new MultiMap<>(new HashMap<>(), () -> new TreeSet<>(SequenceNameLocusComparator.SINGLETON));
    String chr = null;
    for (SequenceNameLocus r : d1.regions().values()) {
      if (chr == null || !chr.equals(r.getSequenceName())) {
        chr = r.getSequenceName();
        chrs.add(chr);
      }
      unionByChr.put(chr, r);
    }
    int chrpos = -1;
    for (SequenceNameLocus r : d2.regions().values()) {
      chr = r.getSequenceName();
      if (chrpos == -1 || !chrs.get(chrpos).equals(chr)) {
        for (int newchrpos = chrpos + 1; newchrpos < chrs.size(); newchrpos++) { // Advance to position of this chromosome if it exists
          if (chrs.get(newchrpos).equals(chr)) {
            chrpos = newchrpos;
            break;
          }
        }
        if (chrpos == -1 || !chrs.get(chrpos).equals(chr)) { // Not found, so need to add
          // But check earlier in the list --> error due to disagreeing chr sorting between the datasets
          for (int newchrpos = chrpos - 1; newchrpos >= 0; newchrpos--) {
            if (chrs.get(newchrpos).equals(chr)) {
              throw new IllegalArgumentException("Cannot UnionJoin datasets: disagreement about ordering of chromosome: " + chr);
            }
          }
          //System.err.println("Inserting new chromosome from dataset: " + chr);
          chrs.add(++chrpos, chr);
        }
      }
      unionByChr.put(chr, r);
    }
    final List<SequenceNameLocus> union = new ArrayList<>();
    chrs.forEach(c -> union.addAll(unionByChr.get(c)));
    return union;
  }

  // Pass through dataset, adding rows as required to ensure regions match
  private void expandRegions(Collection<SequenceNameLocus> regions, RegionDataset d) {
    int i = 0;
    for (SequenceNameLocus s : regions) {
      final SequenceNameLocus r = i < d.size() ? d.regions().get(i) : null;
      if (r == null || SequenceNameLocusComparator.SINGLETON.compare(s, r) != 0) {
        //System.err.println("Inserting row at: " + s);
        d.add(i, s);
      }
      i++;
    }
  }
}
