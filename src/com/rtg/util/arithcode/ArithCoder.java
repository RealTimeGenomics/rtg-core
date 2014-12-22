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

import com.reeltwo.jumble.annotations.TestClass;

/**
 * A generic arithmetic coding subclass containing elements common to
 * both arithmetic decoding and arithmetic coding.
 *
 * @version 1.1
 */
@TestClass("com.rtg.util.arithcode.ArithTest")
class ArithCoder {

  /**
   * The low bound on the current interval for coding.  Initialized
   * to zero.
   */
  protected long mLow; // implied = 0;

  /**
   * The high bound on the current interval for coding.
   * Initialized to top value possible.
   */
  protected long mHigh = TOP_VALUE;

  /**
   * Construct a generic arithmetic coder.
   */
  protected ArithCoder() { }

  /**
   * Precision of coding, expressed in number of bits used for
   * arithmetic before shifting out partial results.
   */
  protected static final int CODE_VALUE_BITS = 27;

  /**
   * The largest possible interval value. All <code>1</code>s.
   */
  protected static final long TOP_VALUE = ((long) 1 << CODE_VALUE_BITS) - 1;

  /**
   * 1/4 of the largest possible value plus one.
   */
  protected static final long FIRST_QUARTER = TOP_VALUE / 4 + 1;

  /**
   * 1/2 of the largest possible value; <code>2 * FIRST_QUARTER</code>
   */
  protected static final long HALF = 2 * FIRST_QUARTER;

  /**
   * 3/4 of the largest possible value; <code>3 * FIRST_QUARTER</code>
   */
  protected static final long THIRD_QUARTER = 3 * FIRST_QUARTER;

}
