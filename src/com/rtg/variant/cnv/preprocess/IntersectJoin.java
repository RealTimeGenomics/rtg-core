/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
 * Joins two datasets dropping non matching regions
 */
public class IntersectJoin extends SimpleJoin {

  IntersectJoin(File in) throws IOException {
    super(in);
  }

  IntersectJoin(RegionDataset dataset, String prefix) {
    super(dataset, prefix);
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    // Drop rows where sequences mismatch
    final Set<String> sequences = getSequenceNames(dataset);
    final Set<String> seq2 = getSequenceNames(mDataset);
    final int initial = sequences.size();
    sequences.retainAll(seq2);
    if (sequences.size() != initial) {
      filterSequences(dataset, sequences);
    }
    if (sequences.size() != seq2.size()) {
      filterSequences(mDataset, sequences);
    }

    // Drop rows where regions mismatch
    filterRegions(mDataset, dataset);

    super.process(dataset);
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
