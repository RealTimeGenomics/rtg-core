package com.rtg.variant.bayes.multisample.cancer;

/**
 * @author Sean A. Irvine
 */
public interface Scorer {
  void add(Double score, int refCount, int altCount);

  int size();

  long getTotalRefCount();

  long getTotalAltCount();
}
