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

import com.rtg.mode.Frame;
import com.rtg.position.output.PositionOutput;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class MockPositionOutput extends IntegralAbstract implements PositionOutput {

  @Override
  public void endAll() {
  }

  @Override
  public void endPosition() {
  }

  @Override
  public void endQuery() {
  }

  @Override
  public void endQuerySequence() {
  }

  @Override
  public void hit(final int seqId, final int posn) {
  }

  @Override
  public boolean integrity() {
    return false;
  }

  @Override
  public void nextQuery(final Frame frame, final int seqId) {
  }

  @Override
  public void nextSequence(int seqId, final int length) {
    // do nothing
  }

  @Override
  public void nextSequence(int seqId, final int length, final int usedLength, final byte[] sequence) {
    // do nothing
  }

  @Override
  public void setPosition(final int position) {
  }

  @Override
  public double score() {
    return 0;
  }

  @Override
  public void toString(final StringBuilder sb) {
  }

  /**
   * Return a clone of this suitable for use in the reverse frame
   * @return the clone
   */
  @Override
  public PositionOutput reverseClone() {
    return this;
  }
}
