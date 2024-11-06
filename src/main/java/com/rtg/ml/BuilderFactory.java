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
