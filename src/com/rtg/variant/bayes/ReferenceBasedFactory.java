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

package com.rtg.variant.bayes;


/**
 * Specifies factories that make objects depending on the current reference code
 * @param <D> type of object that is created
 */
public interface ReferenceBasedFactory<D> {

  /**
   * Make an object given the code for the reference sequence.
   *
   * @param ref the reference code.
   * @return the model.
   */
  D make(int ref);

}
