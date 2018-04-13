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

import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusComparator;

/**
 * Joins two datasets containing exactly matching regions
 */
public class SimpleJoin implements DatasetProcessor {

  protected final String mPrefix;
  protected final RegionDataset mDataset;

  SimpleJoin(File in) throws IOException {
    final String[] split = StringUtils.split(in.toString(), '=', 2);
    final File file = new File(split[0]);
    if (split.length > 1 && split[1].length() > 0) {
      mPrefix = split[1].trim() + "_";
    } else {
      mPrefix = "";
    }
    mDataset = RegionDataset.readFromBed(file);
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   */
  public SimpleJoin(RegionDataset dataset, String prefix) {
    mDataset = dataset;
    mPrefix = prefix;
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    // Check regions are compatible
    final ObjectColumn<SequenceNameLocus> rc1 = dataset.regions();
    final ObjectColumn<SequenceNameLocus> rc2 = mDataset.regions();
    final int count = Math.max(rc1.size(), rc2.size());
    for (int i = 0; i < count; ++i) {
      if (i >= rc1.size() || i >= rc2.size()) {
        throw new NoTalkbackSlimException("Cannot join datasets as the number of rows differ: " + dataset.size() + " != " + mDataset.size());
      }
      final SequenceNameLocus r1 = rc1.get(i);
      final SequenceNameLocus r2 = rc2.get(i);
      if (SequenceNameLocusComparator.SINGLETON.compare(r1, r2) != 0) {
        throw new NoTalkbackSlimException("Cannot join datasets: region differerence at row " + i + "\nA: " + dataset.getBedRecord(i) + "\nB: " + mDataset.getBedRecord(i));
      }
    }

    // Move the columns over
    for (int i = 0; i < mDataset.columns(); ++i) {
      final Column col = mDataset.column(i);
      col.setName(mPrefix + col.getName());
      dataset.addColumn(col);
    }
  }
}
