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
package com.rtg.util;

import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Used when there is an invalid command line parameter.
 */
public class InvalidParamsException extends NoTalkbackSlimException {

  /**
   * Can't be bothered with an ErrorType? Do I have the constructor for YOU!
   * @param message error message
   */
  public InvalidParamsException(String message) {
    super(message);
  }

  /**
   *
   * @param error the error to be raised.
   * @param params parameters for the error message.
   */
  public InvalidParamsException(final ErrorType error, final String... params) {
    super(error, params);
  }
}
