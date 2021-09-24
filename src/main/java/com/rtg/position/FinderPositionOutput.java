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
package com.rtg.position;

import com.rtg.index.Finder;
import com.rtg.position.output.PositionOutput;

/**
 * Has a <code>Finder</code> and a <code>PositionOutput</code>
 */
public class FinderPositionOutput {

  private final Finder mFinder;
  private final PositionOutput mPositionOutput;

  /**
   * Constructs a <code>FinderPositionOutput</code>
   * @param finder the finder
   * @param positionOutput the position output
   */
  public FinderPositionOutput(final Finder finder, final PositionOutput positionOutput) {
    this.mFinder = finder;
    this.mPositionOutput = positionOutput;
  }

  public Finder getFinder() {
    return mFinder;
  }

  public PositionOutput getPositionOutput() {
    return mPositionOutput;
  }


}
