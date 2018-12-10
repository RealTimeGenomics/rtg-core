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
  protected final boolean mFirst;
  protected final boolean mAllowDuplicateColumns;

  SimpleJoin(File in) throws IOException {
    this(RegionDataset.readFromBed(getFile(in)), getPrefix(in));
  }

  private static File getFile(File in) {
    return new File(StringUtils.split(in.toString(), '=', 2)[0]);
  }

  private static String getPrefix(File in) {
    final String[] split = StringUtils.split(in.toString(), '=', 2);
    return (split.length > 1 && split[1].length() > 0) ? split[1].trim() + "_" : "";
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   */
  public SimpleJoin(RegionDataset dataset, String prefix) {
    this(dataset, prefix, false, false);
  }

  /**
   * Constructor
   * @param dataset the dataset to join
   * @param prefix a prefix applied to the joined columns
   * @param addFirst if true, any joined columns appear first
   * @param allowDuplicateColumns if false, do not add any column which already exists in the destination dataset (based on name, incorporating prefix)
   */
  public SimpleJoin(RegionDataset dataset, String prefix, boolean addFirst, boolean allowDuplicateColumns) {
    mDataset = dataset;
    mPrefix = prefix;
    mFirst = addFirst;
    mAllowDuplicateColumns = allowDuplicateColumns;
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    simpleJoin(mDataset, dataset);
  }

  protected void simpleJoin(RegionDataset src, RegionDataset dest) {
    // Check regions are compatible
    final ObjectColumn<SequenceNameLocus> rc1 = dest.regions();
    final ObjectColumn<SequenceNameLocus> rc2 = src.regions();
    final int count = Math.max(rc1.size(), rc2.size());
    for (int i = 0; i < count; ++i) {
      if (i >= rc1.size() || i >= rc2.size()) {
        throw new NoTalkbackSlimException("Cannot join datasets as the number of rows differ: " + dest.size() + " != " + src.size());
      }
      final SequenceNameLocus r1 = rc1.get(i);
      final SequenceNameLocus r2 = rc2.get(i);
      if (SequenceNameLocusComparator.SINGLETON.compare(r1, r2) != 0) {
        throw new NoTalkbackSlimException("Cannot join datasets: region difference at row " + i + "\nA: " + dest.getBedRecord(i) + "\nB: " + src.getBedRecord(i));
      }
    }

    // Move the columns over
    for (int i = 0, j = 0; i < src.columns(); ++i) {
      final Column col = src.column(i);
      final String newName = mPrefix + col.getName();
      if (mAllowDuplicateColumns || dest.columnId(newName) == -1) {
        col.setName(newName);
        if (mFirst) {
          dest.addColumn(j++, col);
        } else {
          dest.addColumn(col);
        }
      }
    }
  }
}
