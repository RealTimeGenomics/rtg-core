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

import com.rtg.launcher.HashingRegion;

/**
 * Pulls together all the hash functionality.
 */
public interface NgsHashFunction extends ReadHashFunction, TemplateHashFunction {


  /**
   * Create a clone of this hash function that is suitable for running in an independent
   * thread. This method will ensure that any shared configuration is copied or that
   * thread-safe access is provided.
   * @param region the sequence region beyond which results should be filtered.
   * @return a clone.
   * @throws IOException If an I/O error occurs
   */
  NgsHashFunction threadClone(HashingRegion region) throws IOException;

  /**
   * Close this thread.
   * @throws IOException If an I/O error occurs
   */
  void threadFinish() throws IOException;

}
