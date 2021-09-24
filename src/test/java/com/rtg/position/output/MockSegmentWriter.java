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
package com.rtg.position.output;




/**
 */

class MockSegmentWriter implements SegmentWriter {

  private final StringBuilder mSB = new StringBuilder();

  @Override
  public void write(final Segment segment, final int searchPosition) {
    mSB.append(segment.seqId()).append("\t").append(segment.start()).append("\t").append(segment.end() - segment.start() + 1).append(com.rtg.util.StringUtils.LS);
  }

  @Override
  public String toString() {
    return mSB.toString();
  }
}
