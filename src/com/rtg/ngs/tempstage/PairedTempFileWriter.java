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
   * @throws IOException if an I/O Error occurs
   */
  void nextTemplateId(long templateId) throws IOException;

  /**
   * Tells the alignment writer where mated results should be written.
   * Must be called before any calls to <code>pairResult</code>
   * @param matedOut output stream to send mated results to.
   * @exception IOException if an I/O error occurs
   */
  void initialiseMated(OutputStream matedOut) throws IOException;

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
   * @exception IOException if an I/O error occurs
   */
  void initialiseUnmated(OutputStream unmatedOut, MapQScoringReadBlocker unmatedBlockerLeft, MapQScoringReadBlocker unmatedBlockerRight) throws IOException;

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
