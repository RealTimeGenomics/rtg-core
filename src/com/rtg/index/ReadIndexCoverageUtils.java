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
package com.rtg.index;

import com.rtg.index.hash.ngs.ReadDecoder;
import com.rtg.util.StringUtils;

/**
 * Find out what reads aren't in the index at all
 */
public final class ReadIndexCoverageUtils {

  private ReadIndexCoverageUtils() {

  }

  /** True if read rejection should be reported */
  //public static final boolean REPORT_REJECTED_READS = Boolean.valueOf(System.getProperty("report.rejected.reads", "false"));

  /** Reinsert rejected reads once identified */
  //public static final boolean REINSERT_REJECTED_READS = Boolean.valueOf(System.getProperty("reinsert.rejected.reads", "false"));

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
