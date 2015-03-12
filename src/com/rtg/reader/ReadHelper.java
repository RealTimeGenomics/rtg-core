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
package com.rtg.reader;

import java.io.IOException;


/**
 * Utility for reading... reads... into a byte array
 */
public final class ReadHelper {

  private ReadHelper() { }

  /**
   * Get a read.
   *
   * @param r reader
   * @param readId read identifier
   * @return read
   */
  public static byte[] getRead(SequencesReader r, long readId) {
    if (r == null) {
      return null;
    }
    try {
      final byte[] b = new byte[r.length(readId)];
      r.read(readId, b);
      return b;
    } catch (final IOException ex) {
      // This should not occur in MemorySequencesReader
      throw new IllegalStateException(ex.getMessage());
    }
  }

  /**
   * Get the quality information for a read.
   *
   * @param r reader
   * @param readId read
   * @return quality
   */
  public static byte[] getQual(SequencesReader r, long readId) {
    try {
      if (r == null || !r.hasQualityData()) {
        return null;
      }
      final byte[] b = new byte[r.length(readId)];
      r.readQuality(readId, b);
      return b;
    } catch (final IOException ex) {
      // This should not occur in MemorySequencesReader
      throw new IllegalStateException(ex.getMessage());
    }
  }

}
