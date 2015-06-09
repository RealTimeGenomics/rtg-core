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
package com.rtg.sam;

import java.util.List;
import java.util.NoSuchElementException;

import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Iterator for skipping records not within target region. Assumes given iterator is sorted.
 */
public class SamRestrictingIterator implements CloseableIterator<SAMRecord> {

  private final ReferenceRanges<String> mRegions;
  private final CloseableIterator<SAMRecord> mIterator;

  private int mTemplate = -1;
  private int mEndTemplate = -1;
  private List<RangeList.RangeData<String>> mSeqRegions = null;
  private int mCurrentRegionIdx = -1;
  private RangeList.RangeData<String> mCurrentRegion = null;

  private SAMRecord mNextRecord;

  /**
   * @param iterator source of records
   * @param regions the restriction regions
   */
  public SamRestrictingIterator(CloseableIterator<SAMRecord> iterator, ReferenceRanges<String> regions) {
    if (regions == null || iterator == null) {
      throw new NullPointerException();
    }
    //System.out.println("SamRestrictingIterator to ranges: " + regions);
    mRegions = regions;
    mIterator = iterator;
    for (int templateId : regions.sequenceIds()) {
      mEndTemplate = templateId;
    }
    populateNext();
  }

  private void populateNext() {
    mNextRecord = null;
    while (mIterator.hasNext()) {
      final SAMRecord rec = mIterator.next();
      final int refId = rec.getReferenceIndex();

      if (refId != mTemplate) { // Get current sequence ranges

        if (refId > mEndTemplate) {  // Past the last template in requested list so we can bail out
          //System.out.println("Aborting iterator for record at " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + " as it is after last template");
          break;
        }

        mTemplate = refId;
        final RangeList<String> seqList = mRegions.get(mTemplate);
        mSeqRegions = (seqList == null || seqList.getRangeList().size() == 0) ? null : seqList.getRangeList();
        mCurrentRegionIdx = 0;
        mCurrentRegion = mSeqRegions == null ? null : mSeqRegions.get(0);
      }

      if (mSeqRegions == null) { // No regions on this template
        continue;
      }

      final int alignmentStart = rec.getAlignmentStart() - 1; // to 0-based
      int alignmentEnd = rec.getAlignmentEnd(); // to 0-based exclusive = noop
      if (alignmentEnd == 0) {
        // Use the read length to get an end point for unmapped reads
        alignmentEnd = rec.getAlignmentStart() + rec.getReadLength();
      }

      // Advance active region if need be
      while (mCurrentRegion.getEnd() <= alignmentStart && mCurrentRegionIdx < mSeqRegions.size() - 1) {
        mCurrentRegionIdx++;
        mCurrentRegion = mSeqRegions.get(mCurrentRegionIdx);
      }

      if (alignmentEnd <= mCurrentRegion.getStart()) {  // before region
        continue;
      }

      if (alignmentStart >= mCurrentRegion.getEnd()) {  // past region
        if (mTemplate == mEndTemplate) {
          //System.out.println("Aborting iterator for record at " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + " as it is after last region");
          break;
        } else {
          continue;
        }
      }

      mNextRecord = rec;
      break;
    }
  }


  @Override
  public void close() {
    mIterator.close();
  }

  @Override
  public boolean hasNext() {
    return mNextRecord != null;
  }

  @Override
  public SAMRecord next() {
    if (mNextRecord == null) {
      throw new NoSuchElementException();
    }
    final SAMRecord ret = mNextRecord;
    populateNext();
    return ret;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

}
