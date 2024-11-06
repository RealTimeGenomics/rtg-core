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
package com.rtg.index.queue;

import com.rtg.index.Add;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.array.longindex.LongChunks;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Provides access to one queue for each radix.
 */
public class IndexQueue extends IntegralAbstract implements Add {

  /*
   * The data in here is represented in a very non-object oriented way because we are
   * being very careful about avoiding cache misses while creating and accessing the queues.
   * The data is referenced by pointers into mMemory from mQueueInfo which has
   * one set of data for each radix value. This set of information consists of:
   * <pre>
   * start      the last position in the first block for the queued items.
   * current    the current position within the block which will be written to next.
   * currentEnd the last position in the current block.
   * length     the total length of the current block.
   * </pre><br>
   * Each block once it has finished being written looks as follows:
   * <pre>
   * if length > 0:
   *    <1 ... length-2> entries
   *    next    the last position in the next block (takes 2 slots)
   *    length  the total length of the block
   * if length < 0:
   *    entries
   *    used    the index of the last used entry in the block (takes 2 slots).
   *    length  the negative of the total length of the block
   * </pre>
   */

  private static final int INTEGER_MAX_BITS = 31;

  /**
   * A power of 2 greater than or equal to the longest cache-line we expect to meet.
   * i7 processors have a cache line of 512 = 8 * 64;
   */
  private static final int MIN_BLOCK_SIZE = 128;

  //Offsets into the individual "fields" for each radix.
  private static final int START =  0;
  private static final int CURR  =  1;
  private static final int END   =  2;
  private static final int LENGTH = 3;
  private static final int TOTAL_BITS  = 2;

  private final int mLowerBits;

  private final int mUpperBits;

  private final int mRadix;

  private final long[] mQueueInfo;

  private final ExtensibleIndex mMemory;

  private boolean mClosed = false;

  private final long mHashMask;

  private final int mValueBits;

  /**
   * @param lowerBits number of bits to right of radix.
   * @param upperBits number of bits in a radix (high order bits of hash).
   * @param length an initial estimate of the number of entries to be stored in the queue.
   *               the memory will be expanded as necessary later, however, initial memory sufficient for this
   *               is allocated.
   * @param valueBits number of value bits.
   */
  public IndexQueue(final int lowerBits, final int upperBits, final long length, final int valueBits) {
    if (upperBits < 0 || lowerBits < 0 || valueBits < 0 || length < 0) {
      throw new IllegalArgumentException();
    }
    final long bs = (length >> upperBits) + 3L;
    if (bs < 2) { // Extremely unlikely to happen, requires length very close to Long.MAX_VALUE and upperBits == 0
      throw new IllegalArgumentException();
    }
    final long blockSize = Math.max(5L, bs); //ensure space for at least two items plus next/length values
    mLowerBits = lowerBits;
    mUpperBits = upperBits;
    mRadix = 1 << mUpperBits;
    mQueueInfo = new long[mRadix << TOTAL_BITS];

    for (int i = 0; i < mRadix; ++i) {
      final int j = i << TOTAL_BITS;
      final long ib = i * blockSize;
      final long st = ib + blockSize - 1;
      mQueueInfo[j + START]  = st;
      mQueueInfo[j + CURR]   = ib;
      mQueueInfo[j + END]    = st;
      mQueueInfo[j + LENGTH] = blockSize;
    }

    //dumpQueueInfo();
    final long memSize = blockSize * mRadix;
    //System.err.println("length=" + length + " upperBits=" + upperBits + " memSize=" + memSize);
    mValueBits = valueBits;
    mHashMask = (1L << mLowerBits) - 1;
    mMemory = makeMemory(memSize);
    assert integrity();
  }

  final ExtensibleIndex makeMemory(final long length) {
    final int lengthBits = MathUtils.ceilPowerOf2Bits(length - 1);
    final int initLengthBits = lengthBits - 7;
    final int initSubarrayBits = initLengthBits <= 0 ? lengthBits : initLengthBits;
    final int subarrayBits;
    if (initSubarrayBits > INTEGER_MAX_BITS) {
      subarrayBits = INTEGER_MAX_BITS;
    } else if (initSubarrayBits <= 8) {
      subarrayBits = 8;
    } else {
      subarrayBits = initSubarrayBits;
    }
    final int bb = lengthBits - subarrayBits;
    final int blockBits = bb < 0 ? 0 : bb;
    return makeMemory(blockBits, length, subarrayBits);
  }

  ExtensibleIndex makeMemory(final int blockBits, final long length, final int subarrayBits) {
    final int maxBits = Math.max(mLowerBits, mValueBits);
    if (maxBits <= 32) {
      return new IntChunks(1 << (blockBits + 1), length, subarrayBits);
    } else {
      return new LongChunks(1 << (blockBits + 1), length, subarrayBits);
    }
  }

  @Override
  public void add(final long hash, final long id) {
    if (mClosed) {
      throw new IllegalStateException();
    }
    final int radix = (int) (hash >>> mLowerBits);
    assert radix < mRadix;
    add(radix, hash & mHashMask);
    add(radix, id);
  }

  private void setFree(final long pos, final long free) {
    if (mMemory instanceof IntChunks) {
      mMemory.setSigned(pos - 1, free >> 32);
      mMemory.set(pos, (int) free);
    } else {
      assert mMemory instanceof LongChunks;
      mMemory.setSigned(pos, free);
    }
  }

  private void add(final int radix, final long value) {
    //System.err.println("radix=" + radix);
    final int i = radix << TOTAL_BITS;
    final long curr = mQueueInfo[i + CURR];
    final long end = mQueueInfo[i + END];
    if (curr < end - 2) { // Last three slots are used for housekeeping
      mMemory.set(curr, value);
      mQueueInfo[i + CURR] = curr + 1;
    } else {
      //System.err.println("Expanding memory block for radix " + radix);
      //need to create a new block
      final long free = mMemory.extendBy(MIN_BLOCK_SIZE);
      final long length = mQueueInfo[i + LENGTH];
      mQueueInfo[i + CURR] = free;
      final long fr = free + MIN_BLOCK_SIZE - 1;
      mQueueInfo[i + END] = fr;
      mQueueInfo[i + LENGTH] = MIN_BLOCK_SIZE;
      setFree(end - 1, fr);
      mMemory.setSigned(end, length);
      // Now that we have made more space, we can be sure that on the following recursion the value will be stored.
      add(radix, value);
    }
  }

  /**
   * Ensure all queues are closed.
   * Should be called before any iterators are returned.
   */
  void close() {
    for (int i = 0; i < mRadix; ++i) {
      final int j = i << TOTAL_BITS;
      final long end = mQueueInfo[j + END];
      final long curr = mQueueInfo[j + CURR];
      //System.err.println("end=" + end + " curr=" + curr);
      setFree(end - 1, curr);
      //mMemory.setSigned(end - 1, curr);
      final long length = mQueueInfo[j + LENGTH];
      mMemory.setSigned(end, -length);
    }
    Diagnostic.developerLog("Index queue memory usage : " + StringUtils.commas(mMemory.bytes()) + " bytes");
    mClosed = true;
  }

  /**
   * Get an iterator over one of the queues.
   * @param radix specifies which queue to return.
   * @param radixShifted radix mask
   * @return the iterator.
   */
  QueueIterator iterator(final int radix, final long radixShifted) {
    if (!mClosed) {
      throw new IllegalStateException();
    }
    final long start = mQueueInfo[(radix << TOTAL_BITS) + START];
    return new QueueIterator(mMemory, start, radixShifted);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IndexQueue lower=").append(mLowerBits).append(" upper=").append(mUpperBits);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final long memLength = mMemory.length();
    for (int i = 0; i < mRadix; ++i) {
      final int j = i << TOTAL_BITS;
      //System.err.println("i=" + i + " start=" + start);
      long st = mQueueInfo[j + START];
      while (true) {
        final long length = mMemory.getSigned(st);
        if (length == 0) {
          break;
        }
        final long us = mMemory.getSigned(st - 1);
        //System.err.println("st=" + st + " length=" + length + " us=" + us);
        if (length < 0) {
          Exam.assertTrue(mClosed);
          Exam.assertTrue("st=" + st + " length=" + length + " us=" + us, st + length + 1 <= us && us < st);
          break;
        } else {
          Exam.assertTrue("length=" + length + " st=" + st + " us=" + us + " memLength=" + memLength, st < us && us < memLength);
          st = us;
        }
      }
      final long curr = mQueueInfo[j + CURR];
      final long end = mQueueInfo[j + END];
      final long qLength = mQueueInfo[j + LENGTH];
      Exam.assertTrue(end - qLength <= curr && curr <= end - 2);
      Exam.assertTrue(0 <= curr && curr < memLength);
      Exam.assertTrue(0 <= end && end < memLength);
      if (mClosed) {
        final long len = mMemory.getSigned(end);
        Exam.assertTrue(len < 0);
        Exam.assertEquals(qLength, -len);
        final long us = mMemory.getSigned(end - 1);
        Exam.assertEquals(us, curr);
      }
    }
    return true;
  }


  @Override
  public final boolean integrity() {
    Exam.assertTrue(mUpperBits + mLowerBits <= 64);
    Exam.assertTrue(mUpperBits < 32);
    Exam.assertEquals(mRadix, 1 << mUpperBits);
    Exam.assertEquals(mRadix << TOTAL_BITS, mQueueInfo.length);
    Exam.assertNotNull(mMemory);
    return true;
  }

}
