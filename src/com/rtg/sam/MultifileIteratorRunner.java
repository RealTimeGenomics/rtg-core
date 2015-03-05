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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.IORunnable;
import com.rtg.util.Populator;
import com.rtg.util.ProgramState;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 * An object to be run in a simple thread pool to get <code>SAMRecords</code>
 * from a set of files, merge them by genome position as it goes, and buffer them
 * for access by the <code>ThreadedMultifileIterator</code>
 * @param <T> type of record to produce
 * Code review done: 28 October 2010 (partial).
 */
@TestClass(value = {"com.rtg.sam.ThreadedMultifileIteratorTest"})
final class MultifileIteratorRunner<T> implements RecordIterator<T>, IORunnable, Comparable<MultifileIteratorRunner<T>> {

  private static final int TIMEOUT = 100000; // milliseconds
  private static final int MAX_BUFFER = 2;  // should be at least two
  private final int mPacketSize; // number of SAM records in each buffer entry

  private final MultifileIterator mIterator;
  private final int mId; // used for tie-breaking
  private final CappedConcurrentLinkedList<Collection<T>> mRecords;
  private final Populator<T> mPopulator;
  private Iterator<T> mPacketIterator = ((Collection<T>) new ArrayList<T>()).iterator();
  private T mTopRecord = null;

  private volatile boolean mVolIsClosing = false;
  private volatile boolean mVolIsFinished = false;

  /**
   * A version of the MultifileIterator to be used with a thread
   * for reading the Sam Records into a queue
   *
   * @param context the SAM reading context
   * @param populator the populator
   * @param id an id number to keep track of instances within a single <code>ThreadedMultifileIterator</code>
   * @param packetSize the size of packets
   * @throws IOException when an error in file access occurs
   */
  MultifileIteratorRunner(SamReadingContext context, Populator<T> populator, int id, int packetSize) throws IOException {
    mIterator = new MultifileIterator(context);
    mId = id;
    mPacketSize = packetSize;
    mPopulator = populator;
    mRecords = new CappedConcurrentLinkedList<>(MAX_BUFFER, mId);
    if (!mIterator.hasNext()) {
      mRecords.setHasNext(false);
    }
  }

  @Override
  public void run() {
    List<T> packet = new ArrayList<>(mPacketSize);
    while (!mVolIsClosing && mIterator.hasNext()) {
      try {
        final SAMRecord next = mIterator.next();
        final boolean hasNext = mIterator.hasNext();
        packet.add(mPopulator.populate(next));
        if (packet.size() >= mPacketSize || !hasNext) {
          mRecords.add(packet, hasNext);
          ProgramState.checkAbort();
          if (hasNext) {
            packet = new ArrayList<>(mPacketSize);
          } else {
            break;
          }
        }
      } catch (final RuntimeException e) {
        synchronized (this) {
          mVolIsFinished = true;
          notifyAll();
        }
        mRecords.setHasNext(false);
        //close();
        throw e;
      }
    }
    synchronized (this) {
      mVolIsFinished = true;
      notifyAll();
    }
  }

  @Override
  public boolean hasNext() {
    return mTopRecord != null || mPacketIterator.hasNext() || !mRecords.isEmpty();
  }

  @Override
  public T next() {
    tryTopNotNull();
    assert mTopRecord != null;
    final T res = mTopRecord;
    if (mPacketIterator.hasNext()) {
      mTopRecord = mPacketIterator.next();
    } else {
      mTopRecord = null;
      tryTopNotNull();
    }
    return res;
  }

  private void tryTopNotNull() {
    if (mTopRecord != null) {
      return;
    }
    if (!mPacketIterator.hasNext()) {
      if (!mRecords.isEmpty()) {
        final Collection<T> packet = mRecords.poll();
        if (packet == null) {
          synchronized (this) {
            mVolIsFinished = true;
            notifyAll();
          }
          throw new IllegalStateException("Top packet in queue is null, the worker threads have probably died");
        }
        mPacketIterator = packet.iterator();
      } else {
        return;
      }
    }
    mTopRecord = mPacketIterator.next();
  }

  @Override
  public int compareTo(MultifileIteratorRunner<T> that) {
    this.tryTopNotNull();
    that.tryTopNotNull();
    final T thisRecord = this.mTopRecord;
    final T thatRecord = that.mTopRecord;

    if (thisRecord == null && thatRecord == null) {
      final int idComp = compareId(that);
      assert idComp != 0 || this == that;
      return idComp;
    }
    if (thisRecord == null) {
      return -1;
    }
    if (thatRecord == null) {
      return +1;
    }

    if (thisRecord instanceof SAMRecord) {
      final int samCompare = SamCompareUtils.compareSamRecords((SAMRecord) thisRecord, (SAMRecord) thatRecord);
      if (samCompare != 0) {
        return samCompare;
      }

    } else if (thisRecord instanceof Comparable) {
      @SuppressWarnings("unchecked")
      final int c = ((Comparable<T>) thisRecord).compareTo(thatRecord);
      if (c != 0) {
        return c;
      }
    }

    final int idCom = compareId(that);
    assert idCom != 0 || this == that;
    return idCom;
  }

  private int compareId(MultifileIteratorRunner<T> that) {
    final int thisId = this.mId;
    final int thatId = that.mId;
    //System.err.println("id  this=" + thisId + " that=" + thatId);
    if (thisId < thatId) {
      return -1;
    }
    if (thisId > thatId) {
      return +1;
    }
    return 0;
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof MultifileIteratorRunner)) {
      return false;
    }
    return mIterator.equals(((MultifileIteratorRunner<?>) obj).mIterator);
  }

  @Override
  public int hashCode() {
    return mIterator.hashCode() + mId;
  }

  @Override
  public synchronized void close() throws IOException {
    mVolIsClosing = true;
    mRecords.close();
    while (!mVolIsFinished) {
      try {
        wait(TIMEOUT);
      } catch (final InterruptedException e) {
        ProgramState.checkAbort();
        throw new IllegalStateException("Interrupted but program not aborting?", e);
      }
    }
    mIterator.close();
  }

  @Override
  public SAMFileHeader header() {
    return mIterator.header();
  }

  @Override
  public void remove() {
    mIterator.remove();
  }

  @Override
  public long getTotalNucleotides() {
    return mIterator.getTotalNucleotides();
  }

  @Override
  public long getInvalidRecordsCount() {
    return mIterator.getInvalidRecordsCount();
  }

  @Override
  public long getFilteredRecordsCount() {
    return mIterator.getFilteredRecordsCount();
  }

  @Override
  public long getOutputRecordsCount() {
    return mIterator.getOutputRecordsCount();
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mIterator.getDuplicateRecordsCount();
  }

  @Override
  public long getTotalRecordsCount() {
    return mIterator.getTotalRecordsCount();
  }

}
