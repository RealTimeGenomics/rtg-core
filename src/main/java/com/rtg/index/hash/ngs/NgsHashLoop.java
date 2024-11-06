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

import java.io.IOException;

import com.rtg.launcher.ISequenceParams;

/**
 * Interface for hash calls
 */
public interface NgsHashLoop {

  /**
   * Scans the buffer and makes calls to
   * <code>hashFunction</code>.
   * @param params specifies reader and start and end positions.
   * @param hashFunction it is given the succesive codes and does the work from there on.
   * @param encoder the read id encoder.
   * @param reverse if true then store the sequence as reverse complement.
   * @return the total number of nucleotides read.
   * @throws IOException if an I/O error occurs.
   */
  long readLoop(final ISequenceParams params, final ReadHashFunction hashFunction, final ReadEncoder encoder, final boolean reverse) throws IOException;

  /**
   * Scans the buffer and makes calls to
   * <code>hashStepTemplate</code> in the hash function.
   * @param params specifies start, end and reader.
   * @param hashFunction does the work as each code is added from the template.
   * @throws IOException if an I/O error occurs.
   */
  void templateLoop(final ISequenceParams params, final TemplateHashFunction hashFunction) throws IOException;

  /**
   * Scans the buffer and makes calls to
   * <code>hashStepTemplate</code> in the hash function.
   * Use threads to run this in parallel
   * @param params specifies start, end and reader.
   * @param hf model used to create on hash function for each thread.
   * @param numberThreads the number of threads to use.
   * @param threadMultiplier thread multiplier
   * @throws IOException if an I/O error occurs.
   */
  void templateLoopMultiCore(final ISequenceParams params, final NgsHashFunction hf, final int numberThreads, int threadMultiplier) throws IOException;
}
