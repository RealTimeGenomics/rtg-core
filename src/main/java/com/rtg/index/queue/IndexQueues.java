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

import java.io.IOException;

import com.rtg.index.Index;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.SizeSplit;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds queues used in multithread index creation.
 */
public class IndexQueues extends IntegralAbstract {

  private static final int DEFAULT_RADIX_BITS = 10;

  private final int mNumberThreads;

  private final IndexQueue[] mQueues;

  private final int mRadixBits;

  private final long mRadixSize;

  private final int mLowerBits;

  /**
   * @param numberThreads number of threads available for parallel execution.
   * @param hashBits number of bits in a hash word.
   * @param size estimated total size of all queues.
   * @param valueBits the number of bits needed to represent a value.
   * @param ipBits number of bits used to generate length of initial pointer table.
   */
  public IndexQueues(final int numberThreads, final int hashBits, final long size, final int valueBits, final int ipBits) {
    if (numberThreads < 1) {
      throw new IllegalArgumentException("threads=" + numberThreads);
    }
    mNumberThreads = numberThreads;
    mRadixSize = (size + numberThreads - 1) / numberThreads;
    mQueues = new IndexQueue[mNumberThreads];
    mRadixBits = radixBits(hashBits, ipBits);
    mLowerBits = hashBits - mRadixBits;
    Diagnostic.developerLog("Lower bits=" + mLowerBits + " radixBits=" + mRadixBits);
    for (int i = 0; i < mQueues.length; ++i) {
      mQueues[i] = new IndexQueue(mLowerBits, mRadixBits, mRadixSize, valueBits);
    }
  }

  final int radixBits(final int hashBits, final int ipBits) {
    return Math.min(hashBits, Math.min(ipBits, DEFAULT_RADIX_BITS));
  }

  /**
   * Get the queue for use by thread q.
   * @param q the thread index.
   * @return the queue.
   */
  public IndexQueue queue(final int q) {
    return mQueues[q];
  }

  private final class FreezeThread implements IORunnable {
    private final Index mIndex;

    private final int mStart;

    private final int mEnd;

    /**
     * freeze thread
     */
    private FreezeThread(final Index index, final int start, final int end) {
      mIndex = index;
      mStart = start;
      mEnd = end;
    }


    @Override
    public void run() {
      final OneShotTimer timer = new OneShotTimer("LR_BS_freeze_add_" + Thread.currentThread().getName().replace(" ", "_"));
      for (int i = mStart; i < mEnd; ++i) {
        for (final IndexQueue mQueue : mQueues) {
          final long radixShifted = ((long) i) << mLowerBits;
          final QueueIterator it = mQueue.iterator(i, radixShifted);
          while (it.next()) {
            mIndex.add(it.hash(), it.id());
          }
        }
      }
      timer.stopLog();
    }
  }

  /**
   * Transfer all values from queues into the index.
   * Does this twice to accord with compressed memory conventions.
   * @param index to be frozen.
   */
  public void freeze(final Index index) {
    final OneShotTimer freezeCloseTimer = new OneShotTimer("LR_BS_freeze_close");
    for (int i = 0; i < mNumberThreads; ++i) {
      mQueues[i].close();
    }
    freezeCloseTimer.stopLog();
    final OneShotTimer freezeTimer1 = new OneShotTimer("LR_BS_freeze1");
    freeze1(index, "IndexFreezePass1");
    freezeTimer1.stopLog();
    final OneShotTimer freezeTimer2 = new OneShotTimer("LR_BS_freeze2");
    freeze1(index, "IndexFreezePass2");
    freezeTimer2.stopLog();
  }

  private void freeze1(final Index index, String label) {
    final SimpleThreadPool pool = new SimpleThreadPool(mNumberThreads, label, true);
    pool.enableBasicProgress(mNumberThreads);
    final SizeSplit ss = new SizeSplit(1 << mRadixBits, mNumberThreads);

    final OneShotTimer freezeAddTimer = new OneShotTimer("LR_BS_freeze_add");
    for (int i = 0; i < mNumberThreads; ++i) {
      pool.execute(new FreezeThread(index, ss.start(i), ss.start(i + 1)));
    }
    try {
      pool.terminate();
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }
    freezeAddTimer.stopLog();
    final OneShotTimer freezeIndexTimer = new OneShotTimer("LR_BS_freeze_index");
    index.freeze();
    freezeIndexTimer.stopLog();
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IndexQueues: threads=").append(mNumberThreads).append(" radixBits=").append(mRadixBits).append(" radixSize=").append(mRadixSize).append(" lowerBits=").append(mLowerBits);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mNumberThreads >= 1);
    Exam.assertTrue(mRadixBits >= 1);
    Exam.assertEquals(mQueues.length, mNumberThreads);
    return true;
  }

}
