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
package com.rtg.util.array.intindex;

import java.io.IOException;
import java.io.ObjectInputStream;

import com.rtg.util.array.IndexType;

/**
 * Contains the only public ways of constructing an IntIndex.
 */
public final class IntCreate {
  private IntCreate() { // private so cant create an instance of this utility class
  }

  /**
   * Create a new IntIndex of the specified length.
   * Chooses an appropriate implementation depending on the length.
   * @param length number of entries in the IntIndex.
   * @return an IntIndex.
   * @exception NegativeArraySizeException if length negative.
   */
  public static IntIndex createIndex(final long length) {
    if (length < 0) {
      throw new NegativeArraySizeException("Negative length=" + length);
    }
    if (length <= IntIndex.MAX_LENGTH) {
      return new IntArray(length);
    } else {
      return new IntChunks(length);
    }
  }

  /**
   * loads an index saved by {@link IntIndex#save(java.io.ObjectOutputStream)}
   * @param stream stream to load from
   * @return the index stored in the stream
   * @throws IOException if an IO error occurs
   */
  public static IntIndex loadIndex(ObjectInputStream stream) throws IOException {
    final IndexType type = IndexType.values()[stream.readInt()];
    switch (type) {
      case ARRAY:
        return IntArray.loadIndex(stream);
      case CHUNKS:
        return IntChunks.loadIndex(stream);
      default:
        throw new IOException("Unrecognized type");
    }
  }
}

