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
