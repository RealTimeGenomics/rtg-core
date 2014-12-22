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
package com.rtg.util.integrity;

/**
 * Used for classes that have integrity constraints.
 *
 */
public interface Integrity {

  /**
   * Check integrity relationships that operate on just the immediate
   * fields of this object. The implementations will normally be a long
   * sequence of assert statements.
   *
   * @return true so that this can be called from assert statements.
   */
  boolean integrity();


  /**
   * Check integrity relationships that operate on this object and
   * objects accessible from it. The implementations will normally be a
   * long sequence of assert statements. Intended for use with complex
   * data structures such as trees that have non-trivial and expensive
   * to check global integrity constraints.
   *
   * @return true so that this can be called from assert statements.
   */
  boolean globalIntegrity();
}

