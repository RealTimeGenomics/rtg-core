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

import java.util.Properties;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Provides mechanism for creating {@link com.rtg.ml.BuildClassifier}'s by name
 */
@TestClass("com.rtg.ml.BaggedClassifierBuilderTest")
public final class BuilderFactory {
  private BuilderFactory() { }

  /**
   * Enumerates known {@link com.rtg.ml.BuildClassifier}'s and provides ability to construct them
   */
  public enum BuilderType {
    /** corresponds to {@link com.rtg.ml.ZeroRBuilder} */
    ZERO_R {
      @Override
      public BuildClassifier create() {
        return new ZeroRBuilder();
      }
    },
    /** corresponds to {@link com.rtg.ml.RandomTreeBuilder} */
    RANDOM_TREE {
      @Override
      public BuildClassifier create() {
        return new RandomTreeBuilder();
      }
    },
    /** corresponds to {@link com.rtg.ml.BaggedClassifierBuilder} */
    BAGGED {
      @Override
      public BuildClassifier create() {
        return new BaggedClassifierBuilder();
      }
    };

    /**
     * Create an empty BuildClassifier of the appropriate type
     * @return the new BuildClassifier
     */
    public abstract BuildClassifier create();
  }

  /**
   * Instantiate a BuildClassifier by name
   * @param builderName the type name
   * @return the classifier
   */
  public static BuildClassifier create(String builderName) {
    return BuilderType.valueOf(builderName).create();
  }

  /**
   * Instantiate a BuildClassifier by name and set configuration parameters
   * @param builderName the type name
   * @param parameters configuration parameters
   * @return the classifier
   */
  public static BuildClassifier create(String builderName, Properties parameters) {
    final BuildClassifier builder = create(builderName);
    builder.setProperties(parameters);
    return builder;
  }
}
