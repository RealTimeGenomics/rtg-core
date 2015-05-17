/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import java.io.Closeable;
import java.io.IOException;

/**
 *
 */
public interface WriterWrapper extends Closeable {

  /**
   * Method to write a specified sequence to the destination.
   * @param seqId the id of the sequence to write.
   * @param dataBuffer an appropriately sized buffer for data.
   * @param qualityBuffer an appropriately sized buffer for quality.
   * @throws IllegalStateException whenever
   * @throws IOException whenever
   */
  void writeSequence(long seqId, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalStateException, IOException;

  @Override
  void close() throws IOException;
}
