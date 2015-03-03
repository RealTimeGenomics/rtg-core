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

package com.rtg.util.arithcode;

/**
 */
interface DetailedModel extends ArithCodeModel {


  /**
   * Returns the total count for the current context.
   * @return Total count for the current context.
   */
  int totalCount();

  /**
   * Returns the symbol whose interval of low and high counts
   * contains the given count.
   * @param count The given count.
   * @return The symbol whose interval contains the given count.
   */
  int pointToSymbol(int count);

  /**
   * Calculates <code>{low count, high count, total count}</code> for
   * the given symbol in the current context.
   * The cumulative counts
   * in the return must be such that <code>0 &lt;= low count &lt; high
   * count &lt;= total count</code>.
   * <P>
   * This method will be called exactly once for each symbol being
   * encoded or decoded, and the calls will be made in the order in
   * which they appear in the original file.
   * @param symbol The next symbol to decode.
   * @param result Array into which to write range.
   */
  void interval(int symbol, int[] result);
}
