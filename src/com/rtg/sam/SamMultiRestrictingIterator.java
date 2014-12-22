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

import java.io.IOException;
import java.util.NoSuchElementException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.tabix.VirtualOffsets;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.SequenceNameLocus;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.RuntimeIOException;

/**
 * Like SamRestrictingIterator but uses VirtualOffsets to perform separate region retrievals.
 */
@TestClass("com.rtg.sam.SamClosedFileReaderTest")
final class SamMultiRestrictingIterator implements CloseableIterator<SAMRecord> {

  private final BlockCompressedInputStream mStream;
  private final VirtualOffsets mOffsets;
  private final SAMFileHeader mHeader;
  private final boolean mIsBam;
  private final String mLabel;

  private int mCurrentOffset = 0;
  private CloseableIterator<SAMRecord> mCurrentIt;
  private SequenceNameLocus mCurrentRegion;
  private int mPreviousAlignmentStart = Integer.MIN_VALUE;
  private long mPreviousVirtualOffsetStart = Long.MIN_VALUE;
  private int mCurrentTemplate;

  private int mDoubleFetched = 0;
  private SAMRecord mNextRecord;

  SamMultiRestrictingIterator(BlockCompressedInputStream stream, VirtualOffsets offsets, SAMFileHeader header, boolean isBam, String label) throws IOException {
    mStream = stream;
    mOffsets = offsets;
    mHeader = header;
    mIsBam = isBam;
    mLabel = label;

    mCurrentIt = null;
    mCurrentOffset = 0;
    mCurrentTemplate = -1;

    // Set up for first region and if it has no data, skip ahead to find one that does
    populateNext(true);
  }


  // Support the ability to "push-back" a record when we get the first one past the bounds of the current region,
  // Since if we end up re-using the underlying iterator, we need to evaluate that record w.r.t the new region.
  private SAMRecord mBuffered = null;
  private void setBuffered(SAMRecord rec) {
    mBuffered = rec;
  }
  private void clearBuffered() {
    mBuffered = null;
  }
  private SAMRecord nextOrBuffered() {
    if (mBuffered == null) {
      return mCurrentIt.next();
    }
    final SAMRecord ret = mBuffered;
    mBuffered = null;
    return ret;
  }


  private void populateNext(boolean force) {
    final SAMRecord previousRecord = mNextRecord;
    mNextRecord = null;
    if (force) {
      advanceSubIterator();
    }
    while (mCurrentIt != null) {
      if (!mCurrentIt.hasNext()) {  // Only happens when stream is exhausted, so effectively just closes things out.
        advanceSubIterator();
      } else {
        final SAMRecord rec = nextOrBuffered();
        final int refId = rec.getReferenceIndex();

        if (refId > mCurrentTemplate) { // current offset has exceeded region and block overlapped next template
          setBuffered(rec);
          advanceSubIterator();
        } else {

          if (refId < mCurrentTemplate) { // Current block may occasionally return records from the previous template if the block overlaps
            continue;
          }

          final int alignmentStart = rec.getAlignmentStart() - 1; // to 0-based
          int alignmentEnd = rec.getAlignmentEnd(); // to 0-based exclusive = noop
          if (alignmentEnd == 0) {
            // Use the read length to get an end point for unmapped reads
            alignmentEnd = rec.getAlignmentStart() + rec.getReadLength();
          }

          if (alignmentEnd <= mCurrentRegion.getStart()) {  // before region
            continue;
          }

          if (alignmentStart <= mPreviousAlignmentStart) { // this record would have been already returned by an earlier region
            //Diagnostic.developerLog("Ignoring record from earlier block at " + rec.getReferenceName() + ":" + rec.getAlignmentStart());
            mDoubleFetched++;
            if (mDoubleFetched % 100000 == 0) {
              Diagnostic.developerLog("Many double-fetched records for source " + mLabel + " noticed at " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + " in region " + mCurrentRegion + " (skipping through to " + mPreviousAlignmentStart + ")");
            }
            continue;
          }

          if (alignmentStart >= mCurrentRegion.getEnd()) {  // past current region, advance the iterator and record the furtherest we got
            if (previousRecord != null && previousRecord.getReferenceIndex() == mCurrentTemplate) {
              mPreviousAlignmentStart = previousRecord.getAlignmentStart();
            } else {
              mPreviousAlignmentStart = Integer.MIN_VALUE;
            }
            setBuffered(rec);
            advanceSubIterator();
            continue;
          }

          mNextRecord = rec;
          break;
        }
      }
    }
  }

  protected void advanceSubIterator() {
    try {
      if (mCurrentOffset < mOffsets.size()) {
        if (mOffsets.start(mCurrentOffset) != mPreviousVirtualOffsetStart) {
          //Diagnostic.developerLog("After region " + mCurrentRegion + ", opening new reader at offset " + VirtualOffsets.offsetToString(mOffsets.start(mCurrentOffset)) + " for region " + mOffsets.region(mCurrentOffset));
          // Don't the current iterator here as it will axe off the underlying block compressed input stream, oddly only for SAM, not BAM
          clearBuffered();
          mStream.seek(mOffsets.start(mCurrentOffset));
          // Warning: Confusing constructors being used here - the true is only used to force a different constructor
          mCurrentIt = mIsBam ? new SAMFileReader(mStream, mHeader, true).iterator() : new SAMFileReader(mStream, mHeader).iterator();

        //} else {
        //  Diagnostic.developerLog("After region " + mCurrentRegion + ", re-using existing reader for region " + mOffsets.region(mCurrentOffset));

        }

        mCurrentRegion = mOffsets.region(mCurrentOffset);
        final int newTemplate = mHeader.getSequenceIndex(mCurrentRegion.getSequenceName());
        if (newTemplate != mCurrentTemplate) {
          mPreviousAlignmentStart = Integer.MIN_VALUE;
        }
        mPreviousVirtualOffsetStart = mOffsets.start(mCurrentOffset);
        mCurrentTemplate = newTemplate;

        mCurrentOffset++;
      } else {
        closeCurrent();
      }
    } catch (IOException e) {
      throw new RuntimeIOException(e.getMessage(), e);
    }
  }

  @Override
  public void close() {
    closeCurrent();
    Diagnostic.developerLog("There were " + mDoubleFetched + " SAM records double-fetched due to overlapping blocks");
  }

  private void closeCurrent() {
    if (mCurrentIt != null) {
      final CloseableIterator<SAMRecord> currentIt = mCurrentIt;
      mCurrentIt = null;
      currentIt.close();
    }
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
    populateNext(false);
    return ret;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
