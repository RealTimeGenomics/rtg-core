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

import java.util.Arrays;

/**
 * Encapsulate an instance containing both attribute data and binary classification.
 *
 */
public class Instance {

  private final double[] mInstance;
  private final boolean mPositive;
  private double mWeight;

  /**
   * Create the instance
   * @param instance the attribute data
   * @param isPositive true if the instance is a positive example
   */
  public Instance(double[] instance, boolean isPositive) {
    this(instance, isPositive, 1.0);
  }

  /**
   * Create the instance
   * @param instance the attribute data
   * @param isPositive true if the instance is a positive example
   * @param weight the weight of the instance
   */
  public Instance(double[] instance, boolean isPositive, double weight) {
    mInstance = instance;
    mPositive = isPositive;
    mWeight = weight;
  }

  /**
   * Used to handle missing values during training.
   * @param weight weight to assign the copy
   * @return a version of the instance with a new weight
   */
  public Instance reweight(final double weight) {
    return new Instance(mInstance, mPositive, weight);
  }

  /**
   * Override the existing weight, not a copy.
   * @param weight new weight
   */
  void setWeight(final double weight) {
    mWeight = weight;
  }

  /**
   * @return a copy of the instance (attribute values can be altered in this new copy without affecting the original)
   */
  public Instance deepCopy() {
    return new Instance(Arrays.copyOf(mInstance, mInstance.length), mPositive, mWeight);
  }

  /**
   * @return the attribute data
   */
  public double[] instance() {
    return mInstance;
  }

  /**
   * @return true if the instance is a positive example
   */
  public boolean isPositive() {
    return mPositive;
  }

  /**
   * @return the weight of the instance
   */
  public double weight() {
    return mWeight;
  }
}
