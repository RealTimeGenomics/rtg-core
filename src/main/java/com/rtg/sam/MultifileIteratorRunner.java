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
package com.rtg.sam;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.IORunnable;
import com.rtg.util.Populator;
import com.rtg.util.ProgramState;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.WarningType;

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

  private static final class Packet<T> extends ArrayList<T> {
    private boolean mHasNext = false;
    private Packet(int size) {
      super(size);
    }
  }

  private static final int TIMEOUT = 2;
  private static final int MAX_BUFFER = 2;  // should be at least two
  private final int mPacketSize; // number of SAM records in each buffer entry

  private final MultifileIterator mIterator;
  private final int mId; // used for tie-breaking
  private final LinkedBlockingQueue<Packet<T>> mRecords;
  private final Populator<T> mPopulator;
  private Iterator<T> mPacketIterator = new ArrayList<T>().iterator();
  private boolean mMorePackets = true;
  private T mTopRecord = null;

  private long mInvalidRecords;

  private volatile boolean mClosed = false;

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
    mRecords = new LinkedBlockingQueue<>(MAX_BUFFER);
    if (!mIterator.hasNext()) {
      mMorePackets = false;
    }
  }

  @Override
  public void run() {
    Packet<T> packet = new Packet<>(mPacketSize);
    while (!mClosed && mIterator.hasNext()) {
      try {
        final SAMRecord next = mIterator.next();
        final boolean hasNext = mIterator.hasNext();
        final T populate = mPopulator.populate(next);
        if (populate != null) {
          packet.add(populate);
        } else {
          ++mInvalidRecords;
          maybeWarn(next);
        }
        if (packet.size() >= mPacketSize || !hasNext) {
          packet.mHasNext = hasNext;
          mRecords.put(packet);
          ProgramState.checkAbort();
          if (hasNext) {
            packet = new Packet<>(mPacketSize);
          } else {
            break;
          }
        }
      } catch (InterruptedException e1) {
        mClosed = true;
        ProgramState.checkAbort();
        throw new IllegalStateException("Interrupted but program not aborting?", e1);
      } catch (final Throwable t) {
        mClosed = true;
        throw t;
      }
    }
  }

  private void maybeWarn(SAMRecord record) {
    if (mInvalidRecords <= 5) {
      Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING1, record.toString());
      Diagnostic.userLog("Invalid record: " + record);
    }
  }

  @Override
  public boolean hasNext() {
    tryTopNotNull();
    return mTopRecord != null;
  }

  @Override
  public T next() {
    assert mTopRecord != null;
    final T res = mTopRecord;
    mTopRecord = null;
    return res;
  }

  private void tryTopNotNull() {
    if (mTopRecord != null) {
      return;
    }
    if (mPacketIterator.hasNext()) {
      mTopRecord = mPacketIterator.next();
      return;
    }
    while (!mClosed && (mMorePackets || !mRecords.isEmpty())) {
      try {
        final Packet<T> packet = mRecords.poll(TIMEOUT, TimeUnit.SECONDS);
        if (packet != null) {
          mPacketIterator = packet.iterator();
          mMorePackets = packet.mHasNext;
          mTopRecord = mPacketIterator.next();
          return;
        }
      } catch (InterruptedException e) {
        ProgramState.checkAbort();
        throw new IllegalStateException("Interrupted but program not aborting?", e);
      }
    }
    if (mClosed) {
      ProgramState.checkAbort();
      throw new IllegalStateException("Worker thread has aborted");
    }
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
    } else if (thisRecord instanceof SamReaderRecord) {
      final int samCompare = SamCompareUtils.compareSamRecords(((SamReaderRecord) thisRecord).inner(), ((SamReaderRecord) thatRecord).inner());
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
    return Integer.compare(this.mId, that.mId);
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
  public void close() throws IOException {
    mClosed = true;
    mRecords.clear(); // Ensures that if we close before all data is consumed, any pending put() is unblocked
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
    return mIterator.getInvalidRecordsCount() + mInvalidRecords;
  }

  @Override
  public long getFilteredRecordsCount() {
    return mIterator.getFilteredRecordsCount();
  }

  @Override
  public long getOutputRecordsCount() {
    return mIterator.getOutputRecordsCount() - mInvalidRecords;
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mIterator.getDuplicateRecordsCount();
  }

  @Override
  public long getOverCoverageRecordsCount() {
    return mIterator.getOverCoverageRecordsCount();
  }

  @Override
  public long getTotalRecordsCount() {
    return mIterator.getTotalRecordsCount();
  }

}
