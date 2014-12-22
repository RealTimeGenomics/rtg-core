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
package com.rtg.index.hash.ngs;

import java.io.Closeable;
import java.io.IOException;

import com.rtg.launcher.HashingRegion;

/**
 * Abstract class for accumulating and writing outputs.
 *
 */
public interface OutputProcessor extends Closeable {

  /**
   * Process the given results.
   * @param templateName name of template
   * @param frame frame of result
   * @param readId read number (0 based).
   * @param tStart template start (0 based).
   * @param score simple score
   * @param scoreIndel score indel
   * @throws IOException if problems writing output
   */
  void process(long templateName, String frame, int readId, int tStart, int score, int scoreIndel) throws IOException;

  /**
   * Informs the output processor that hit processing is finished and it is time to do any
   * end of run processing/cleanup/summarizing
   * @throws IOException if problems when writing output.
   */
  void finish() throws IOException;

  /**
   * Create a clone of this OutputProcessor that is suitable for running in an independent
   * thread. This method will ensure that any shared configuration is copied or that
   * thread-safe access is provided.
   *
   * @param region a clipping region, to which thread result output must be restricted.
   * @return a clone.
   * @throws IOException If an I/O error occurs
   */
  OutputProcessor threadClone(HashingRegion region) throws IOException;


  /**
   * Finish up any thread local business.
   * @throws IOException if an I/O Error occurs
   */
  void threadFinish() throws IOException;
}

