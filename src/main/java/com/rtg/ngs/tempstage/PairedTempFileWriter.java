/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.ngs.tempstage;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.pairedend.MatedHitInfo;

/**
 * Defines a device for writing pairings to the temp files
 *
 */
public interface PairedTempFileWriter extends Closeable {

  /**
   * Step to the specified template identifier.  The ids should be monotonically
   * increasing. <code>Long.MAX_VALUE</code> should be used to indicate no more sequences.
   *
   * @param templateId a <code>long</code> value
   * @throws IOException if an I/O error occurs.
   */
  void nextTemplateId(long templateId) throws IOException;

  /**
   * Tells the alignment writer where mated results should be written.
   * Must be called before any calls to <code>pairResult</code>
   * @param matedOut output stream to send mated results to.
   */
  void initialiseMated(OutputStream matedOut);

  /**
   * Write a mating for the left (lesser template position) hit of a
   * pair. This method should populate the alignment fields of the hit
   * info object.
   *
   * @param hitInfo a <code>MatedHitInfo</code> containing the hit coordinates.
   * @exception IOException if an error occurs.
   * @return true if the pair was considered a good pairing (the
   * caller can use this to determine whether to call <code>pairResultRight</code>
   * for this pairing)
   */
  boolean pairResultLeft(MatedHitInfo hitInfo) throws IOException;

  /**
   * Write a mating for the right (lesser template position) hit of a
   * pair.
   *
   * @param hitInfo a <code>MatedHitInfo</code> containing the hit coordinates.
   * @exception IOException if an error occurs.
   */
  void pairResultRight(MatedHitInfo hitInfo) throws IOException;

  /**
   * Closes mated output
   * @throws IOException if an IO exception occurs
   */
  void closeMated() throws IOException;


  /**
   * Tells the alignment writer where unmated results should be written.
   * Must be called before any calls to <code>unmatedOut</code>
   * @param unmatedOut output stream to send unmated results to.
   * @param unmatedBlockerLeft blocker to keep track of best score per read on left
   * @param unmatedBlockerRight blocker to keep track of best score per read on right
   */
  void initialiseUnmated(OutputStream unmatedOut, MapQScoringReadBlocker unmatedBlockerLeft, MapQScoringReadBlocker unmatedBlockerRight);

  /**
   * Write an unmated result
   * @param readId read identifier
   * @param first true if read is first in pair
   * @param rc true if match was on reverse strand
   * @param start start position of match
   * @return true if result was written, false if it was filtered
   * @throws IOException if an io exception occurs
   */
  boolean unmatedResult(int readId, boolean first, boolean rc, int start) throws IOException;

  /**
   * Closes unmated output
   * @throws IOException if an IO exception occurs
   */
  void closeUnmated() throws IOException;

}
