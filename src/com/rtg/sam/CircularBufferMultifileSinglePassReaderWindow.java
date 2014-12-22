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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import com.rtg.util.Populator;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Reader window that uses a {@link com.rtg.util.io.ClosedFileInputStream} to keep files closed during normal operation and just works from a buffer.
 * May return overflow records to indicate that records have been discarded.
 * @param <T> record type
 */
public class CircularBufferMultifileSinglePassReaderWindow<T extends ReaderRecord<T> & MateInfo> implements ReaderWindow<T> {

  private static final int DEFAULT_BUFFER_LENGTH = 2;
  private ReaderRecord<T>[] mBuffer; // Primary circular buffer
  private int[] mDepth; // Chain length at this position
  private int[] mMaxLengths; // Max record length at this position
  private int[] mMinimumStart; // Stores offsets for how far back in buffer to search for overlapping records

  // Holds the one record that is ahead of the requested region.
  private T mSingleRecordStash = null;

  private final int mDepthLimit; // Maximum depth of records to retain at any start position
  private final Populator<T> mPopulator;
  private int mFirstPosition; // The leftmost start position coordinate that we may be asked about (we may store some records with start positions left of that) (0 based in template co-ordinates)
  private int mLastPosition; // The rightmost start position coordinate that we may be asked about (?exclusive) (0 based in template co-ordinates)
  private int mFirstStart; // The left most start position of any read (0 based in template co-ordinates, inclusive)
  private int mLastEnd; // The rightmost end position of any read (0 based in template co-ordinates, exclusive) (always >= mFlushedTo)
  private int mMaxSeenDepth = 0;
  private long mDroppedRecords = 0; // Total count of records discarded due to overflow

  private final RecordIterator<T> mIterator;
  private final int mTemplateIndex;
  private List<FlushLocus> mFlushLocus; //union of all disjoint regions that flush has been called on
  private int mFlushedTo; // Fully flushed to this position (both bins and records)

  /**
   * Utility method to create a RecordIterator, for testing purposes.
   *
   * @param <V> the record type
   * @param samFiles files to read records from
   * @param filterParams parameters for filtering
   * @param numReadingThreads no of threads used for reading SAM files
   * @param pop the record populator
   * @return a RecordIterator
   * @throws IOException if an IO error occurs
   */
  public static <V extends ReaderRecord<V> & MateInfo> RecordIterator<V> defaultIterator(Collection<File> samFiles, SamFilterParams filterParams, int numReadingThreads, Populator<V> pop) throws IOException {
    final SingletonPopulatorFactory<V> pf = new SingletonPopulatorFactory<>(pop);
    RecordIterator<V> it = new ThreadedMultifileIterator<>(samFiles, numReadingThreads, pf, filterParams, SamUtils.getUberHeader(samFiles));
    if (filterParams.findAndRemoveDuplicates()) {
      it = new DedupifyingRecordIterator<>(it);
    }
    return it;
  }

  /**
   * Note: the supplied record iterator must be closed explicitly by the caller when finished using it as the
   * close method for the <code>CircularBufferMultifileSinglePassReaderWindow</code> will not close it.
   * @param recordIt iterator supplying records
   * @param pop the record populator
   * @param templateIndex sequence index for template
   * @param templateStart first position being retrieved (zero based)
   * @param maxDepth maximum records per position
   */
  public CircularBufferMultifileSinglePassReaderWindow(RecordIterator<T> recordIt, Populator<T> pop, int templateIndex, int templateStart, int maxDepth) {
    mIterator = recordIt;
    mPopulator = pop;
    mTemplateIndex = templateIndex;
    mFirstPosition = templateStart == -1 ? 0 : templateStart;
    mBuffer = createArray(DEFAULT_BUFFER_LENGTH);
    mDepth = new int[DEFAULT_BUFFER_LENGTH];
    mMaxLengths = new int[DEFAULT_BUFFER_LENGTH];
    mDepthLimit = Math.max(1, maxDepth); // blows up deeper with NPE if 0 is allowed through
    mMinimumStart = new int[DEFAULT_BUFFER_LENGTH];
    mLastPosition = mFirstPosition;
    mFlushedTo = mFirstPosition;
    mFirstStart = mFirstPosition;
    mLastEnd = mFlushedTo;
    mFlushLocus = new ArrayList<>();
    mFlushLocus.add(new FlushLocus(mFirstPosition, mFirstPosition));
  }

  /**
   * Returns the position where data has actually been flush to.
   * @return the position flushed to.
   */
  public int flushedTo() {
    return mFlushedTo;
  }

  /**
   * Returns the start of the earliest chuck that has not had a flush request.
   * @return position of minimum chunk with work outstanding.
   */
  public int finishedTo() {
    return mFlushLocus.get(0).mEnd;
  }

  private void addFlushLocus(int start, int end) {
    mFlushLocus = addFlushLocus(mFlushLocus, new FlushLocus(start, end));
    //dumpFlushList();
  }

  static ArrayList<FlushLocus> addFlushLocus(final List<FlushLocus> flushLocus, final FlushLocus lnew0) {
    final ArrayList<FlushLocus> newlist = new ArrayList<>();
    FlushLocus lnew = lnew0;
    for (final FlushLocus fl : flushLocus) {
      if (fl.isJoinable(lnew)) {
        lnew.join(fl);
      } else if (fl.mStart < lnew.mStart) {
        newlist.add(fl);
      } else {
        newlist.add(lnew);
        lnew = fl;
      }
    }
    newlist.add(lnew);
    return newlist;
  }

  @Override
  public void flush(final int start, final int end) throws IOException {
    //Diagnostic.developerLog("flush start=" + start + " " + end);
    if (start < mFirstPosition) {
      throw new IllegalArgumentException("Flush start " + start + " is less than " + mFirstPosition);
    }
    if (end > mLastPosition) {
      throw new IllegalArgumentException("Flush end " + end + " is greater than " + mLastPosition);
    }

    // Now look at where it's safe to actually flush to
    addFlushLocus(start, end);

    final FlushLocus reallyFlush = mFlushLocus.get(0);
    //Diagnostic.developerLog("reallyFlush: " + reallyFlush);
    final int newStart = mFlushedTo;
    final int overlapEnd = reallyFlush.mEnd; // Will clear overlap offsets up to this position
    //final int binEnd = overlapEnd + mMinimumStart[overlapEnd % mBuffer.length]; // Will clear record buffer up to this position

    final int binEnd = overlapEnd + ((mLastEnd < overlapEnd) ? 0 : mMinimumStart[overlapEnd % mBuffer.length]);
    //System.err.println("binEnd:" + binEnd + " : " + overlapEnd + " : " + mLastEnd + " : " + mBuffer.length + " : " + (overlapEnd % mBuffer.length) + " : " + mMinimumStart[overlapEnd % mBuffer.length]);

    assert binEnd <= overlapEnd;
    if (binEnd <= newStart) {
      return;
    }
    // Clear out records
    if (mSingleRecordStash != null && mSingleRecordStash.getStart() < binEnd) {
      mSingleRecordStash = null;
    }
    for (int i = mFirstStart; i < binEnd; i++) {
      final int ix = i % mBuffer.length;
      mBuffer[ix] = null;
      mDepth[ix] = 0;
      mMaxLengths[ix] = 0;
    }
    for (int i = newStart; i < overlapEnd; i++) {
      final int ix = i % mBuffer.length;
      mMinimumStart[ix] = 0;
    }

    //Diagnostic.developerLog("MinStart.length: " + mMinimumStart.length);
    //Diagnostic.developerLog("mBuffer.length: " + mBuffer.length);
    //Diagnostic.developerLog("overlapEnd: " + overlapEnd);
    //Diagnostic.developerLog("binEnd: " + binEnd);

    mFlushedTo = binEnd;
    mFirstStart = binEnd;
    mFirstPosition = overlapEnd;
    if (mFirstPosition > mLastEnd) {
      mLastEnd = mFirstPosition;
    }
    //Diagnostic.developerLog("CBMSPRW - Discarding start=" + start + " end=" + end + " count=" + count);
    //System.err.println("flush(" + start + "," + end + ") fully flushed to " + mFlushedTo + ", range now " + mFirstPosition + "," + mLastPosition);
  }

  @Override
  public void advanceBuffer(int end) {
    mLastPosition = Math.max(end, mLastPosition);
    if (mSingleRecordStash != null) {
      if (addRecord(mSingleRecordStash)) {
        mSingleRecordStash = null;
      } else {
        return;
      }
    }
    final RecordIterator<T> skippy = mIterator;
    while (skippy.hasNext()) {
      final T r = skippy.next();
      if (!addRecord(r)) {
        mSingleRecordStash = r;
        break;
      }
    }
  }

  //true means continue looping, false means stash record and stop
  private boolean addRecord(T r) {
    final int alignmentStart = r.getStart();
    final int alignmentEnd = alignmentStart + r.getLength();
    if (r.getSequenceId() < mTemplateIndex || (r.getSequenceId() == mTemplateIndex && alignmentEnd <= mFirstPosition)) {
      //not overlapping and before region.
      return true;
    }
    if (alignmentStart >= mLastPosition || r.getSequenceId() > mTemplateIndex) {
      //beyond region
      return false;
    }
    resizeToPosition(alignmentStart, alignmentEnd); // Alignment end to ensure any existing overlap offset doesn't get clobbered
    //System.err.println("add record start=" + alignmentStart + " end=" + alignmentEnd);
    insertRecord(r);
    return true;
  }

  protected void resizeToPosition(final int start, final int end) {
    //Diagnostic.developerLog("resizeToPosition: " + start + " : " + end);
    final int newStart = Math.min(start, mFirstStart);
    final int newEnd = Math.max(end, mLastEnd);
    while (newEnd > newStart + mBuffer.length) {
      resize();
    }
    mFirstStart = newStart;
    mLastEnd = newEnd;
    //Diagnostic.developerLog("resizedTo: " + mBuffer.length);
  }

  private void resize() {
    final int resize = (3 * mBuffer.length + 2) / 2;
    //System.err.println("resize before=" + mBuffer.length + " after=" + resize);
    final ReaderRecord<T>[] newBuffer = createArray(resize);
    final int[] newDepth = new int[resize];
    final int[] newMaxLenths = new int[resize];
    final int[] newMinimumStart = new int[resize];
    //copy the buffers
    for (int i = mFirstStart; i < mLastEnd; i++) {
      final int j = i % mBuffer.length;
      final int k = i % resize;
      newBuffer[k] = mBuffer[j];
      newDepth[k] = mDepth[j];
      newMaxLenths[k] = mMaxLengths[j];
    }

    for (int i = mFirstPosition; i < mLastEnd; i++) {
      final int j = i % mBuffer.length;
      final int k = i % resize;
      newMinimumStart[k] = mMinimumStart[j];
    }
    mBuffer = newBuffer;
    mDepth = newDepth;
    mMaxLengths = newMaxLenths;
    mMinimumStart = newMinimumStart;
  }

  public long getValidRecordsCount() {
    return mIterator.getOutputRecordsCount();
  }

  /**
   * @return number of invalid SAM records read so far
   */
  public long getInvalidRecordsCount() {
    return mIterator.getInvalidRecordsCount();
  }

  public long getFilteredRecordsCount() {
    return mIterator.getFilteredRecordsCount();
  }

  private void insertRecord(T record) {
    final int start = record.getStart();
    final int index = start % mBuffer.length;
    //System.err.println("add buffer start=" + start + " index=" + index);
    final int length = record.getLength();
    mMaxLengths[index] = Math.max(mMaxLengths[index], length);
    if (++mDepth[index] >= mDepthLimit) {
      if (mDepth[index] > mMaxSeenDepth) {
        mMaxSeenDepth = mDepth[index];
      }
      if (mDepth[index] == mDepthLimit) {
        Diagnostic.developerLog("CBMSPRW tossing all records with start position: " + start);
        mBuffer[index] = mPopulator.overflow(start, mMaxLengths[index]);
        mDroppedRecords += mDepthLimit;
      } else {
        mDroppedRecords++;
        if (length > mBuffer[index].getLength()) {
          // Update discard region to be maximum length observed
          mBuffer[index] = mPopulator.overflow(start, length);
        }
      }
    } else {
      @SuppressWarnings("unchecked")
      final T item = (T) mBuffer[index];
      record.setNextInChain(item);
      mBuffer[index] = record;
    }

    // Update overlap information
    for (int delta = 0, i = start; delta < length; delta++, i++) {
      if (i >= mFirstPosition) {
        final int ix = i % mBuffer.length;
        mMinimumStart[ix] = Math.min(-delta, mMinimumStart[ix]);
      }
    }
  }

  //It is worth noting that the definition of the flush method works at odds with this method, given it tells us what we're done with regarding start positions
  //and says nothing about what we're done with regarding overlapping positions
  @Override
  public Iterator<T> recordsOverlap(final int start, final int end) throws IOException {
    //System.err.println("recordsOverlap start=" + start + " end=" + end);
    // Can ask for start == end in the case of looking for records with indels just to the right? of start
    final int correctedEnd0 = start == end ? end + 1 : end;
    if (start < mFirstPosition) {
      throw new IllegalArgumentException("Expected next call to start at >= " + mFirstPosition + " got: " + start);
    }
    if (correctedEnd0 <= start) {
      throw new IllegalArgumentException("End " + end + " must be greater than start " + start + ".");
    }
    advanceBuffer(correctedEnd0);
    final int correctedEnd = Math.min(mLastEnd, correctedEnd0);
    assert end <= mLastPosition;
    final List<T> records = new ArrayList<>();
    final int startIx0 = start % mBuffer.length;
    final int correctedStart = start + mMinimumStart[startIx0]; // Adjust by offset for overlap
    assert correctedStart >= mFirstStart;
    for (int i = correctedStart; i < correctedEnd; i++) {
      final int ix = i % mBuffer.length;
      @SuppressWarnings("unchecked")
      T rec = (T) mBuffer[ix];
      while (rec != null) {
        if (rec.getStart() + rec.getLength() > start) { // Does it overlap the start?
          records.add(rec);
        }
        rec = rec.chain();
      }
    }
    //System.err.println("recordsOverlap size=" + records.size());
    //System.err.println("CBMSPR |overlap(" + start + ", " + end + ")| = " + records.size();
    return records.iterator();
  }

  /**
   * Total number of records discarded due to extreme coverage regions.
   * @return discarded record count
   */
  public long getTossedRecordCount() {
    return mDroppedRecords;
  }

  /**
   * Close underlying inputs.
   */
  public void close() {
    Diagnostic.developerLog("CBMSPRW total dropped records = " + mDroppedRecords);
    Diagnostic.developerLog("CBMSPRW maximum observed overdepth = " + mMaxSeenDepth);
    //mIterator.close();
  }

  private ReaderRecord<T>[] createArray(final int size) {
    @SuppressWarnings("unchecked")
    final ReaderRecord<T>[] ts = (ReaderRecord<T>[]) new ReaderRecord<?>[size];
    return ts;
  }
}
