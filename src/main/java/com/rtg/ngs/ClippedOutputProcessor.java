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
package com.rtg.ngs;

import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;

/**
 * A class that can be constructed via <code>threadClone</code> to provide results
 * clipping. The clipping is based on raw hit positions (rather than
 * post-alignment), so should not be used for OutputProcessors that
 * perform alignment themselves.
 *
 */
public class ClippedOutputProcessor implements OutputProcessor {

  private final OutputProcessor mInner;

  private final HashingRegion mRegion;


  /**
   * Creates a new <code>ClippedOutputProcessor</code> instance.
   *
   * @param inner an <code>OutputProcessor</code> to which non-filtered hits are delegated
   * @param region a <code>ClipRegion</code> value
   */
  public ClippedOutputProcessor(OutputProcessor inner, HashingRegion region) {
    mInner = inner;
    mRegion = region;
  }

  // Implementation of com.rtg.index.hash.ngs.OutputProcessor

  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
    if (mRegion.isInRange(templateId, tStart)) {
      mInner.process(templateId, frame, readId, tStart, score, scoreIndel);
    }
  }

  @Override
  public void finish() {
    // Noop
  }

  @Override
  public OutputProcessor threadClone(final HashingRegion clipRegion) {
    throw new UnsupportedOperationException();
  }
  @Override
  public void threadFinish() throws IOException {
    mInner.threadFinish();
  }

  @Override
  public void close() {
    // Noop
  }

}
