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
package com.rtg.ml;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Provides mechanism for loading {@link PredictClassifier}'s
 */
public final class MlPredictLoader {
  private MlPredictLoader() { }

  /**
   * Enumerates known {@link PredictClassifier}'s and provides ability to load them
   */
  public enum MlPredictType {
    /** corresponds to {@link com.rtg.ml.ZeroRBuilder.ZeroRClassifier} */
    ZERO_R {
      @Override
      public PredictClassifier loadClassifier(DataInputStream dis, Dataset data) throws IOException {
        return new ZeroRBuilder.ZeroRClassifier(dis);
      }
    },
    /** corresponds to {@link BinaryTreeClassifier} */
    BINARY_TREE {
      @Override
      public PredictClassifier loadClassifier(DataInputStream dis, Dataset data) throws IOException {
        return new BinaryTreeClassifier(dis, data);
      }
    },
    /** corresponds to {@link BaggedClassifier} */
    BAGGED {
      @Override
      public PredictClassifier loadClassifier(DataInputStream dis, Dataset data) throws IOException {
        return new BaggedClassifier(dis, data);
      }
    };

    /**
     * loads a classifier indicated by enum value
     * @param dis the stream to load from
     * @param data set of attributes for encoding/decoding values
     * @return the loaded classifier
     * @throws IOException if an IO error occurs
     */
    public abstract PredictClassifier loadClassifier(DataInputStream dis, Dataset data) throws IOException;
  }

  /**
   * Loads the next classifier in the stream
   * @param is stream to load from
   * @param data set of attributes for encoding/decoding values
   * @return the classifier
   * @throws IOException if an IO error occurs
   */
  public static PredictClassifier loadPredictClassifier(InputStream is, Dataset data) throws IOException {
    return loadPredictClassifier(new DataInputStream(is), data);
  }

  /**
   * Loads the next classifier in the stream
   * @param dis stream to load from
   * @param data set of attributes for encoding/decoding values
   * @return the classifier
   * @throws IOException if an IO error occurs
   */
  public static PredictClassifier loadPredictClassifier(DataInputStream dis, Dataset data) throws IOException {
    final int predictType = dis.readInt();
    if (predictType < 0 || predictType >= MlPredictType.values().length) {
      throw new IOException("Prediction type out of range. Could the model be corrupt?");
    }
    final MlPredictType type = MlPredictType.values()[predictType];
    return type.loadClassifier(dis, data);
  }
}
