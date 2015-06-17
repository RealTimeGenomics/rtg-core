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
package com.rtg.variant.bayes.multisample.cancer;

/**
 * Keeps statistical counts on counts of reference and alternate alleles.
 */
public interface Scorer {

  /**
   * Contribute counts.
   * @param score score for record
   * @param refCount reference allele count
   * @param altCount alternate allele count
   */
  void add(Double score, int refCount, int altCount);

  /**
   * Number of items retained.
   * @return total number of records retained
   */
  int size();

  /**
   * @return total count of reference allele for retained records
   */
  long getTotalRefCount();

  /**
   * @return total count of alternate allele for retained records
   */
  long getTotalAltCount();
}
