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

package com.rtg.variant.realign;

/**
 * Implementations that use delta algorithm to speed up all-paths calculations.
 */
public interface DeltaInterface extends AllPaths {

  /**
   * Specify the hypothesis to replace in the template while leaving
   * read and rest of template and initial all-paths calculations undisturbed.
   * @param replace string to be replaced in template.
   */
  void setReplace(String replace);

}
