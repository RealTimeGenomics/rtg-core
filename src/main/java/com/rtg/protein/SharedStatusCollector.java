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
package com.rtg.protein;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.mode.TranslatedFrame;
import com.rtg.ngs.MapStatistics;
import com.rtg.reader.Arm;
import com.rtg.ngs.MapStatisticsField;
import com.rtg.reader.NamesInterface;
import com.rtg.util.StringUtils;
import com.rtg.util.io.IOUtils;

/**
 * This class stores the status of individual reads,
 * and caches the protein-and-frame version of each read.
 *
 */
public class SharedStatusCollector {
  //Sync locks for status array in superclass object
  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  private static final byte TAB = (byte) '\t';
  private static final byte[] LS = StringUtils.LS.getBytes();

  static final byte EXCEEDS_ALIGNMENT_THRESHOLD = 0x01;
  static final byte EXCEEDS_PERCENT_ID_THRESHOLD = 0x02;
  static final byte EXCEEDS_E_SCORE_THRESHOLD = 0x04;
  static final byte BELOW_BIT_SCORE_THRESHOLD = 0x08;
  static final byte EXCEEDS_N_THRESHOLD = 0x10;
  static final byte RESULT_WRITTEN = 0x20;
  static final byte FAILED_THRESHOLD_FLAGS = 0x1F;

  //cast below is for c-sharp
  private static final byte EXCEEDS_ALIGNMENT_THRESHOLD_CHAR = (byte) 'd';
  private static final byte EXCEEDS_N_THRESHOLD_CHAR = (byte) 'e';
  private static final byte EXCEEDS_PERCENT_ID_THRESHOLD_CHAR = (byte) 'f';
  private static final byte EXCEEDS_E_SCORE_THRESHOLD_CHAR = (byte) 'g';
  private static final byte EXCEEDS_BIT_SCORE_THRESHOLD_CHAR = (byte) 'h';

  // stores the status of the reads
  private final byte[] mReadsStatus;

  /* Array of locks for multiple threads */
  private final Object[] mThreadLocks;

  /**
   * A cache for the protein version of each read-frame combination.
   * This is indexed by read id and frame number.
   */
  private final byte[][] mReadCacheProtein;

  private final MapStatistics mStatistics;

  SharedStatusCollector(int numberOfReads, MapStatistics statistics) {
    mReadsStatus = new byte[numberOfReads];
    mStatistics = statistics;
    mReadCacheProtein = new byte[numberOfReads * TranslatedFrame.values().length][];
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  void setStatus(int readId, byte status) {
    synchronized (mThreadLocks[readId & THREAD_LOCK_MASK]) {
      mReadsStatus[readId] |= status;
    }
  }

  /**
   * @param r read id and frame number
   * @return read translated into protein space, or null if not translated yet.
   */
  protected final byte[] getReadProtein(int r) {
    return mReadCacheProtein[r];
  }

  /**
   * Initialize the protein version of a read for a given frame.
   * This can be used by multiple threads, provided they all put the
   * same data in there, and it remains constant after it has been set.
   *
   * @param r read id and frame number
   * @param protein the read translated into protein.  Must be non-null.
   */
  protected final void putReadProtein(int r, byte[] protein) {
    mReadCacheProtein[r] = protein;
  }

  void writeUnmapped(final OutputStream unmapped, NamesInterface readsNames, long readIdOffset) throws IOException {
    for (int i = 0; i < mReadsStatus.length; ++i) {
      final byte status = mReadsStatus[i];
      if ((status & RESULT_WRITTEN) == RESULT_WRITTEN) {
        continue;
      }
      if (readsNames == null) {
        IOUtils.writeLong(unmapped, readIdOffset + i); //we know its not written so write it in unmapped
      } else {
        unmapped.write(readsNames.name(i).getBytes());
      }
      if ((status & EXCEEDS_N_THRESHOLD) == EXCEEDS_N_THRESHOLD) {
        unmapped.write(TAB);
        unmapped.write(EXCEEDS_N_THRESHOLD_CHAR);
      } else if ((status & BELOW_BIT_SCORE_THRESHOLD) == BELOW_BIT_SCORE_THRESHOLD) {
        unmapped.write(TAB);
        unmapped.write(EXCEEDS_BIT_SCORE_THRESHOLD_CHAR);
      } else if ((status & EXCEEDS_E_SCORE_THRESHOLD) == EXCEEDS_E_SCORE_THRESHOLD) {
        unmapped.write(TAB);
        unmapped.write(EXCEEDS_E_SCORE_THRESHOLD_CHAR);
      } else if ((status & EXCEEDS_PERCENT_ID_THRESHOLD) == EXCEEDS_PERCENT_ID_THRESHOLD) {
        unmapped.write(TAB);
        unmapped.write(EXCEEDS_PERCENT_ID_THRESHOLD_CHAR);
      } else if ((status & EXCEEDS_ALIGNMENT_THRESHOLD) == EXCEEDS_ALIGNMENT_THRESHOLD) {
        unmapped.write(TAB);
        unmapped.write(EXCEEDS_ALIGNMENT_THRESHOLD_CHAR);
      }
      unmapped.write(LS);
    }
  }

  protected void calculateStatistics() {
    if (mStatistics != null) {
      for (final byte status : mReadsStatus) {
        mStatistics.increment(MapStatisticsField.TOTAL_READS, Arm.LEFT);
        if ((status & RESULT_WRITTEN) != 0) {
          mStatistics.increment(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT);
          continue;
        }
        if ((status & FAILED_THRESHOLD_FLAGS) != 0) {
          mStatistics.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT);
        } else {
          mStatistics.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
        }
      }
    }
  }
  protected byte getStatus(int read) {
    return mReadsStatus[read];
  }

}
