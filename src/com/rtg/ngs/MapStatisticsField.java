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
package com.rtg.ngs;

/**
*/
public enum MapStatisticsField {
  /** Key for statistic total number of reads, expected value in form {@link Long} */
  TOTAL_READS,
  /** Key for statistic total number of unique mated reads, expected value in form {@link Long} */
  MATED_UNIQUE_READS,
  /** Key for statistic total number of ambiguous mated reads, expected value in form {@link Long} */
  MATED_AMBIG_READS,
  /** Key for statistic total number of uniquely unmated but mapped reads, expected value in form {@link Long} */
  UNMATED_UNIQUE_READS,
  /** Key for statistic total number of ambiguously unmated but mapped reads, expected value in form {@link Long} */
  UNMATED_AMBIG_READS,

  //unmapped stats
  /** Key for statistic total number of unmapped reads, expected value in form {@link Long} */
  UNMAPPED_NO_HITS,
  /** Key for statistic total number of unmapped reads due to blocking, expected value in form {@link Long} */
  UNMAPPED_BLOCKED,
  /** Key for statistic total number of unmapped reads due to xc d, expected value in form {@link Long} */
  UNMAPPED_MATED_POOR,
  /** Key for statistic total number of unmapped reads due to xc e, expected value in form {@link Long} */
  UNMAPPED_MATED_TOO_MANY,
  /** Key for statistic total number of unmapped reads due to filtering in topN (xc C), expected value in form {@link Long} */
  UNMAPPED_TOPN,
  /** Key for statistic total number of unmapped reads due to xc D, expected value in form {@link Long} */
  UNMAPPED_UNMATED_POOR,
  /** Key for statistic total number of unmapped reads due to xc E, expected value in form {@link Long} */
  UNMAPPED_UNMATED_TOO_MANY,

  //Bad stats
  /** Key for statistic total number of missing reads from first in pair, expected value in form {@link Long} */
  MISSING
}
