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
