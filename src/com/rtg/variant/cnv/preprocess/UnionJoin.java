/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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
    RegionDataset src = mDataset.copy();
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
