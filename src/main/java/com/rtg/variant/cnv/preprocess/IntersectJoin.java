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
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.IntervalComparator;
import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Joins two datasets dropping non matching regions.
 */
public class IntersectJoin extends SimpleJoin {

  IntersectJoin(File in) throws IOException {
    super(in);
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   */
  public IntersectJoin(RegionDataset dataset, String prefix) {
    super(dataset, prefix);
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   * @param addFirst if true, any joined columns appear first
   * @param allowDuplicateColumns if false, do not add any column which already exists in the destination dataset (based on name, incorporating prefix)
   */
  public IntersectJoin(RegionDataset dataset, String prefix, boolean addFirst, boolean allowDuplicateColumns) {
    super(dataset, prefix, addFirst, allowDuplicateColumns);
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    // Drop rows where sequences mismatch
    final RegionDataset src = mDataset.copy();
    final Set<String> sequences = getSequenceNames(dataset);
    final Set<String> seq2 = getSequenceNames(src);
    final int initial = sequences.size();
    sequences.retainAll(seq2);
    if (sequences.size() != initial) {
      filterSequences(dataset, sequences);
    }
    if (sequences.size() != seq2.size()) {
      filterSequences(src, sequences);
    }

    // Drop rows where regions mismatch
    filterRegions(src, dataset);

    simpleJoin(src, dataset);
  }

  // Filter dataset to only keep rows with regions in common
  private void filterRegions(RegionDataset d1, RegionDataset d2) {
    // Step through both sets removing non-matching rows
    int i1 = 0;
    int i2 = 0;
    final ObjectColumn<SequenceNameLocus> rc1 = d1.regions();
    final ObjectColumn<SequenceNameLocus> rc2 = d2.regions();
    SequenceNameLocus last1 = null;
    SequenceNameLocus last2 = null;
    while (i1 < rc1.size() && i2 < rc2.size()) {
      final SequenceNameLocus r1 = rc1.get(i1);
      final SequenceNameLocus r2 = rc2.get(i2);
      if (!r1.getSequenceName().equals(r2.getSequenceName())) {
        if (r1.getSequenceName().equals(last1.getSequenceName())) {
          d1.remove(i1);
          continue;
        } else {
          d2.remove(i2);
          continue;
        }
      }
      if (last1 != null && last1.getSequenceName().equals(r1.getSequenceName()) && IntervalComparator.SINGLETON.compare(r1, last1) < 0) {
        throw new NoTalkbackSlimException("Dataset is not sorted: " + r1.toString() + " < " + last1.toString());
      }
      last1 = r1;
      if (last2 != null && last2.getSequenceName().equals(r2.getSequenceName()) && IntervalComparator.SINGLETON.compare(r2, last2) < 0) {
        throw new NoTalkbackSlimException("Dataset is not sorted: " + r2 + " < " + last2);
      }
      last2 = r2;

      final int c = IntervalComparator.SINGLETON.compare(r1, r2);
      if (c < 0) {
        d1.remove(i1);
      } else if (c > 0) {
        d2.remove(i2);
      } else {
        i1++;
        i2++;
      }
    }
    for (int i = d1.size() - 1; i >= i1; --i) {
      d1.remove(i);
    }
    for (int i = d2.size() - 1; i >= i2; --i) {
      d2.remove(i);
    }
  }

  // Filter dataset to only keep rows against the supplied sequence names
  private void filterSequences(RegionDataset dataset, Set<String> sequences) {
    final ObjectColumn<SequenceNameLocus> rc1 = dataset.regions();
    for (int i = dataset.size() - 1; i >= 0; i--) {
      if (!sequences.contains(rc1.get(i).getSequenceName())) {
        dataset.remove(i);
      }
    }
  }

  private Set<String> getSequenceNames(RegionDataset dataset) {
    final Set<String> seq1 = new HashSet<>();
    dataset.regions().values().forEach(r -> seq1.add(r.getSequenceName()));
    return seq1;
  }
}
