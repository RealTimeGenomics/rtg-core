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

import java.util.HashMap;
import java.util.Iterator;

import com.rtg.util.Utils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Removes duplicate reads from underlying iterator.
 * A read is considered a duplicate if it has the same start position and mate start position as another read.
 * If it is single end it is considered a duplicate if it has the same start and end position as another read.
 * @param <T> record type
 */
public class DedupifyingIterator<T extends ReaderRecord<T> & MateInfo> implements Iterator<T> {

  private static class ReadDuplicateAttributes<X extends ReaderRecord<X> & MateInfo> implements Comparable<ReadDuplicateAttributes<?>> {
    private final boolean mMated;
    private final X mRec;

    ReadDuplicateAttributes(X rec) {
      mRec = rec;
      mMated = rec.isMated();
    }

    private static int compareInts(int a, int b) {
      return Integer.compare(a, b);
    }

    int firstRef() {
      return mRec.getSequenceId();
    }

    int firstPos() {
      return mRec.getStart();
    }

    int secondRef() {
      return mMated ? mRec.getMateSequenceId() : mRec.getSequenceId();
    }

    int secondPos() {
      return mMated ? mRec.getStart() + mRec.getFragmentLength() : mRec.getStart() + mRec.getLength();
    }

    public X record() {
      return mRec;
    }

    /**
     * Records coming from the file should be strictly increasing in this sort order.
     * @param o item to compare to
     * @return like a regular compare method
     */
    public int incomingSortOrderCompareTo(ReadDuplicateAttributes<?> o) {
      final int ret = compareInts(firstRef(), o.firstRef());
      if (ret != 0) {
        return ret;
      }
      return compareInts(firstPos(), o.firstPos());
    }

    //Duplicates must be considered equal according to this comparator method
    @Override
    public int compareTo(ReadDuplicateAttributes<?> o) {
      final int gen = compareInts(mRec.getGenome(), o.mRec.getGenome());
      if (gen != 0) {
        return gen;
      }
      final int ret = incomingSortOrderCompareTo(o);
      if (ret != 0) {
        return ret;
      }
      if (mMated != o.mMated) {
        return mMated ? -1 : 1;
      }
      if (mMated) {
        final int mret = compareInts(secondRef(), o.secondRef());
        if (mret != 0) {
          return mret;
        }
        return compareInts(secondPos(), o.secondPos());
      } else {
        return System.identityHashCode(this) - System.identityHashCode(o); //unmated can never be equal
      }
    }

    @Override
    public boolean equals(Object o) {
      return o instanceof ReadDuplicateAttributes && compareTo((ReadDuplicateAttributes<?>) o) == 0;
    }

    @Override
    public int hashCode() {
      return Utils.pairHash(firstRef(), firstPos(), secondRef(), secondPos(), mMated ? 1 : 0);
    }

  }

  private final Iterator<T> mWrapped;
  private ReadDuplicateAttributes<T> mNext;
  private final DedupifyingList<ReadDuplicateAttributes<T>> mItemsNew;

  // Used when the reference sequence changes to stash the first alignment from the next reference sequence
  private ReadDuplicateAttributes<T> mOverrun;
  protected long mNumDeduped;


  /**
   * @param it iterator to wrap
   */
  public DedupifyingIterator(Iterator<T> it) {
    mWrapped = it;
    mItemsNew = new DedupifyingList<>();
    mNumDeduped = 0;
    populateNext();
  }

  private void addRecord(final ReadDuplicateAttributes<T> rec) {
    // There is no guarantee from the underlying iterator that records will be presented
    // in exactly the same order, so we have to take care here to break ties in a
    // consistent manner. See Bug#1431.
    if (rec.mMated) {
      if (mItemsNew.contains(rec)) {
        final ReadDuplicateAttributes<T> r = mItemsNew.get(rec);
        if (rec.mRec.disambiguateDuplicate(r.mRec) < 0) {
          mItemsNew.replace(rec);
        }
        ++mNumDeduped;
      } else {
        mItemsNew.add(rec, true);
      }
    } else {
      //since unmated is never deduped is wasteful to put in hash set
      //we noticed when running with high coverage single end the run slowed to a crawl previously
      mItemsNew.add(rec, false);
    }
  }

  private void fillItems() {
    assert mItemsNew.isEmpty();
    ReadDuplicateAttributes<T> first = mOverrun;
    if (mOverrun != null) {
      addRecord(mOverrun);
      mOverrun = null;
    }
    while (mWrapped.hasNext()) {
      final T rec = mWrapped.next();
      final ReadDuplicateAttributes<T> current = new ReadDuplicateAttributes<>(rec);
      if (first == null) {
        first = current;
      }
      final int cmp = first.incomingSortOrderCompareTo(current);
      if (cmp == 0) {
        addRecord(current);
      } else if (cmp < 0) {
        mOverrun = current;
        break;
      } else if (current.firstRef() == -1) {
        // Always treat unplaced unmapped as overrun.
        // This means they just get passed through immediately in the order that they are appear.
        // (We don't need to give them a deterministic ordering since they are not being deduplicated.)
        mOverrun = current;
        break;
      } else {
        throw new NoTalkbackSlimException("SAM file was not in sort order:\n" + first.mRec.toString() + "\n" + current.mRec.toString());
      }
    }
  }

  private void populateNext() {
    if (mItemsNew.isEmpty()) {
      fillItems();
    }
    mNext = null;
    final ReadDuplicateAttributes<T> first = mItemsNew.removeFirst();
    if (first != null) {
      mNext = first;
    }
  }

  @Override
  public boolean hasNext() {
    return mNext != null;
  }

  @Override
  public T next() {
    try {
      return mNext.record();
    } finally {
      populateNext();
    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported.");
  }

  /**
   * Put stuff in, get back out in order. We remove duplicates if asked on insertion
   * Not actually a java list
   * @param <T> type of ingoing stuff
   */
  static class DedupifyingList<T> {

    private DedupWrapper<T> mFirst;
    private DedupWrapper<T> mLast;
    private final HashMap<T, DedupWrapper<T>> mMap;

    DedupifyingList() {
      mMap = new HashMap<>();
    }

    /**
     * @param t item to add
     * @param dedup item is to be added to duplicate removal list
     */
    public void add(T t, boolean dedup) {
      assert t != null;
      final DedupWrapper<T> wrap = new DedupWrapper<>(t);
      if (dedup) {
        if (mMap.put(t, wrap) != null) {
          //this is here since the consistent tie-breaking has to be handled externally
          throw new IllegalArgumentException("Don't call this method with duplicates, call replace instead");
        }
      }
      addList(wrap);
    }

    /**
     * replace equivalent value
     * @param t replacement value
     */
    public void replace(T t) {
      //deduping is implied
      assert t != null;
      final DedupWrapper<T> wrap = new DedupWrapper<>(t);
      final DedupWrapper<T> old = mMap.put(t, wrap);
      assert old != null;
      if (old.mPrev != null) {
        old.mPrev.mNext = old.mNext;
      }
      if (old.mNext != null) {
        old.mNext.mPrev = old.mPrev;
      }
      if (mFirst == old) {
        mFirst = old.mNext;
      }
      if (mLast == old) {
        mLast = old.mPrev;
      }
      addList(wrap);
    }

    /**
     * @return true if list is empty
     */
    public boolean isEmpty() {
      return mFirst == null;
    }

    /**
     * @param t key value
     * @return true if <code>t</code> is in the duplicate removal set
     */
    public boolean contains(T t) {
      return mMap.containsKey(t);
    }

    /**
     * Get a stored value from the duplicate removal set
     * @param t key for the value
     * @return value that equals <code>t</code> stored in the list
     */
    public T get(T t) {
      return mMap.get(t).mVal;
    }

    /**
     * Return the first item on the list and remove it
     * @return the value at the start of the list
     */
    public T removeFirst() {
      final DedupWrapper<T> ret = mFirst;
      if (ret == null) {
        return null;
      }
      mMap.remove(ret.mVal);
      mFirst = ret.mNext;
      if (ret.mNext != null) {
        ret.mNext.mPrev = null;
      }
      if (mLast == ret) {
        mLast = null;
      }
      return ret.mVal;
    }

    /**
     * Add to end of list
     * @param wrap item to add
     */
    private void addList(DedupWrapper<T> wrap) {
      if (mFirst == null) {
        mFirst = wrap;
      }
      if (mLast != null) {
        mLast.mNext = wrap;
        wrap.mPrev = mLast;
      }
      mLast = wrap;
    }
  }

  private static class DedupWrapper<T> {
    private final T mVal;
    private DedupWrapper<T> mPrev;
    private DedupWrapper<T> mNext;
    DedupWrapper(T val) {
      mVal = val;
    }
  }

  /**
   * @return number of duplicate records removed
   */
  public long numFiltered() {
    return mNumDeduped;
  }
}
