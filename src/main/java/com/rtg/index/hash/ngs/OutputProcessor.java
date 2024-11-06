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

