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

import java.io.IOException;

/**
 * Handles writing information from a <code>Segment</code> when it is ready to be
 * flushed.
 */
public interface SegmentWriter {

  /**
   * Write the information for this segment.
   * @param segment to be written.
   * @param searchPosition search position (start position of last window of hit)
   * @throws IOException if error while doing write.
   */
  void write(final Segment segment, int searchPosition) throws IOException;

}
