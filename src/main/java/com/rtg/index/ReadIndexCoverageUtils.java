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
package com.rtg.index;

import com.rtg.index.hash.ngs.ReadDecoder;
import com.rtg.util.StringUtils;

/**
 * Find out what reads aren't in the index at all
 */
public final class ReadIndexCoverageUtils {

  private ReadIndexCoverageUtils() {
  }

  /**
   * Find reads which are not present in any of the indexes
   * @param index array of indexes to search
   * @param numReads total number of reads/arms, (specifically max sequence id present in index + 1)
   * @param seqIdDivider what the internal index value should be divided by to acquire the arm id
   * @return array containing id's of arms that are not present in the indexes
   */
  public static int[] findRejectedReads(final Index[] index, final int numReads, final int seqIdDivider) {
    final byte[] readsStatus = new byte[numReads];
    for (final Index anIndex : index) {
      for (long i = 0; i < anIndex.numberEntries(); ++i) {
        final long value = anIndex.getValue(i);
        final int seqId = (int) (value / seqIdDivider);
        readsStatus[seqId] = 1;
      }
    }
    int tot = 0;
    for (final byte readsStatu : readsStatus) {
      tot += readsStatu == 0 ? 1 : 0;
    }
    final int[] ret = new int[tot];
    int k = 0;
    for (int i = 0; i < readsStatus.length; ++i) {
      if (readsStatus[i] == 0) {
        ret[k++] = i;
      }
    }
    return ret;
  }

  /**
   * Summarise the list of rejected arm id's
   * @param rejectedReads as from <code>findRejectedReads</code>
   * @param decode handle to convert arm id to read id
   * @return the summary
   */
  public static String summariseRejectedReads(final int[] rejectedReads, final ReadDecoder decode) {
    final StringBuilder sb = new StringBuilder();
    sb.append("Rejected ").append(decode == ReadDecoder.PAIRED_END ? "arms" : "reads").append(" total: ").append(rejectedReads.length).append(StringUtils.LS);
    for (final int rejectedRead : rejectedReads) {
      sb.append("Rejected: ").append(decode.decode(rejectedRead)).append(" first: ").append(decode.isFirst(rejectedRead)).append(StringUtils.LS);
    }
    return sb.toString();
  }

}
