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
package com.rtg.simulation.reads;

import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 */
public class FragmentTooSmallException extends NoTalkbackSlimException {

  /**
   * @param fragmentLength The length of the fragment the read was to be generated from
   * @param readLength the target read length
   */
  public FragmentTooSmallException(int fragmentLength, int readLength) {
    super(String.format("Fragment length (%d) was not large enough to accomodate read of length %d after indels, try increasing fragment length", fragmentLength, readLength));
  }
}
