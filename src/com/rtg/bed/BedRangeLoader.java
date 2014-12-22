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

package com.rtg.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RangeList.RangeData;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;

/**
 */
@TestClass("com.rtg.vcf.VcfAnnotatorCliTest")
public abstract class BedRangeLoader {

  private final ReferenceRanges.Accumulator mRangeData = new ReferenceRanges.Accumulator();

  private final int mMinAnnotations;
  private int mExtendWarningCount = 0;

  /**
   *
   * @param minAnnotations the minimum number of annotation columns each record in the bed file must contain
   */
  protected BedRangeLoader(int minAnnotations) {
    mMinAnnotations = minAnnotations;
  }

  /**
   * Loads annotations from a set of bed files.
   * @param bedFiles the bed files containing regions to load
   * @throws IOException if an exception occurs while reading a bed file
   */
  public void loadRanges(Collection<File> bedFiles) throws IOException {
    for (final File bedFile : bedFiles) {
      loadRanges(bedFile);
    }
    if (mExtendWarningCount >= 10) {
      Diagnostic.warning("Zero length range extension occurred " + mExtendWarningCount + " times.");
    }
  }

  /**
   * Loads annotations from a bed file.
   * @param bedFile the bed file containing regions to load
   * @throws IOException if an exception occurs while reading a bed file
   */
  public void loadRanges(File bedFile) throws IOException {
    try (BufferedReader in = new BufferedReader(new InputStreamReader(FileUtils.createInputStream(bedFile, true)))) {
      try (BedReader reader = new BedReader(in, mMinAnnotations)) {
        loadRanges(reader);
      }
    }
  }

  private void loadRanges(BedReader reader) throws IOException {
    while (reader.hasNext()) {
      final BedRecord rec = reader.next();
      mRangeData.addRangeData(rec.getSequenceName(), getRangeData(rec));
    }
  }

  /**
   * Merges overlaps between regions and returns the result as a map of <code>RangeList</code> objects, keyed by chromosome name.
   *
   * @return a map of <code>RangeList</code> objects keyed by chromosome name.
   */
  public ReferenceRanges getReferenceRanges() {
    return mRangeData.getReferenceRanges();

  }

  /**
   * Get a range data item corresponding to the given BED record.
   * @param rec the bed record.
   * @return the range data item, or null if this record should be skipped
   */
  protected RangeData<String> getRangeData(BedRecord rec) {
    final int start = rec.getStart();
    int end = rec.getEnd();
    if (end == start) {
      end++;
      // warning - need to have a range of at least 1
      if (mExtendWarningCount < 10) {
        Diagnostic.warning("Zero length range, extending end by 1 : " + rec.toString());
      }
      mExtendWarningCount++;
    }

    return new RangeData<>(start, end, getMeta(rec));
  }

  /**
   * Get the <code>metadata</code> annotation for the specified bed record
   * @param rec the bed record to return a <code>metadata</code> annotation from
   * @return a <code>metadata</code> annotation
   */
  protected abstract String getMeta(BedRecord rec);

}
