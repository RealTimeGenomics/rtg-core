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
package com.rtg.util.format;

/**
 */
public interface FormattedValue {

  /**
   * Describe in a string buffer.
   *
   * @param sb string buffer.
   */
  void toString(StringBuilder sb);


  /**
   * The maximum length any item in this class will return in its
   * <code>toString()</code> method. If there is no known upper bound then
   * Integer.MAX_VALUE is returned.
   *
   * @return maximum length any item in this class will return
   */
  int maxLength();

}

