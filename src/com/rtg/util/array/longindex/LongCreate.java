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
package com.rtg.util.array.longindex;

import java.io.IOException;
import java.io.ObjectInputStream;

import com.rtg.util.array.IndexType;

/**
 * Contains the only public ways of constructing a LongIndex.
 */
public final class LongCreate {
  private LongCreate() { // private so cant create an instance of this utility class
  }

  /**
   * Create a new LongIndex of the specified length.
   * Chooses an appropriate implementation depending on the length.
   * @param length number of entries in the LongIndex.
   * @return a LongIndex.
   * @exception NegativeArraySizeException if length negative.
   */
  public static LongIndex createIndex(final long length) {
    if (length < 0) {
      throw new NegativeArraySizeException("Negative length=" + length);
    }
    // SAI: It seems it is not always possible to get exactly Integer.MAX_VALUE
    // array entries.  Perhaps the JVM uses some slots for housekeeping.
    if (length <= Integer.MAX_VALUE - 5) {
      return new LongArray(length);
    } else {
      return new LongChunks(length);
    }
  }

  /**
   * Create extensible long array
   * @return the array
   */
  public static LongChunks createExtensibleIndex() {
    return new LongChunks(0, 20); //8MiB per chunk
  }

  /**
   * loads an index saved by {@link LongIndex#save(java.io.ObjectOutputStream)}
   * @param stream stream to load from
   * @return the index stored in the stream
   * @throws IOException if an IO error occurs
   */
  public static LongIndex loadIndex(ObjectInputStream stream) throws IOException {
    final IndexType type = IndexType.values()[stream.readInt()];
    switch (type) {
      case ARRAY:
        return LongArray.loadIndex(stream);
      case CHUNKS:
        return LongChunks.loadIndex(stream);
      default:
        throw new IOException("Unrecognized type");
    }
  }
}

