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
 * Stores hits to be processed later. It is up to the implementation which results are kept and which are disposed of.
 */
public interface UptoNStore {

  /**
   * Potentially add hit to list.
   * @param templateId sequence id of reference sequence
   * @param reverse true if hit is reverse complement
   * @param encodedReadId identifier for read
   * @param tStart position on reference of hit
   * @param scoreIndel score for hit
   */
  void process(long templateId, boolean reverse, int encodedReadId, int tStart, int scoreIndel);

  /**
   * Add hits for given read id to results
   * @param results results container
   * @param encodedReadId identifier for read
   */
  void setResults(MatchResult results, int encodedReadId);

  /**
   * @return Debug information about data structure usage
   */
  String histogram();
}
