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
  MISSING,

  /** Key for statistic total number of ignored(short) reads, expected value in form {@link Long} */
  IGNORED,

}
