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
