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

import java.util.ArrayList;

import com.rtg.util.PortableRandom;

/**
 * Encapsulate a dataset suitable for training a two-class classifier.
 *
 */
public class Dataset {

  private final Attribute[] mAttr;
  private double mPosWeight = 0;
  private double mNegWeight = 0;
  private long mPosCount = 0;
  private final ArrayList<Instance> mInstances = new ArrayList<>();

  /**
   * Construct a container for a dataset with specified attributes.
   * @param attr attributes
   */
  public Dataset(final Attribute... attr) {
    assert attr.length > 0;
    mAttr = attr;
  }

  /**
   * @return a copy of the dataset (instances can be added or altered without affecting the original)
   */
  public Dataset deepCopy() {
    final Dataset result = new Dataset(mAttr);
    result.mNegWeight = mNegWeight;
    result.mPosWeight = mPosWeight;
    result.mInstances.ensureCapacity(mInstances.size());
    for (Instance inst : mInstances) {
      result.mInstances.add(inst.deepCopy());
    }
    return result;
  }

  /**
   * Add a training instance. It is assumed the fields of the instance
   * match those of the attributes.
   * @param instance the instance
   */
  public void addInstance(final Instance instance) {
    assert instance.instance().length == mAttr.length;
    mInstances.add(instance);
    if (instance.isPositive()) {
      mPosWeight += instance.weight();
      ++mPosCount;
    } else {
      mNegWeight += instance.weight();
    }
  }

  /**
   * Add all the training instances from another dataset. The dataset must share the exact same attributes.
   * @param dataset the dataset to add
   */
  public void addDataset(Dataset dataset) {
    if (mAttr != dataset.mAttr) {
      throw new IllegalArgumentException("Can only add datasets that share the same attributes");
    }
    mPosCount += dataset.mPosCount;
    mPosWeight += dataset.mPosWeight;
    mNegWeight += dataset.mNegWeight;
    mInstances.ensureCapacity(mInstances.size() + dataset.size());
    mInstances.addAll(dataset.mInstances);
  }

  public Attribute[] getAttributes() {
    return mAttr;
  }

  public ArrayList<Instance> getInstances() {
    return mInstances;
  }

  /** @return the number of instances in the dataset */
  public int size() {
    return mInstances.size();
  }

  /** @return the total weight of positive instances in the dataset */
  public double totalPositiveWeight() {
    return mPosWeight;
  }

  /** @return the total weight of negative instances in the dataset */
  public double totalNegativeWeight() {
    return mNegWeight;
  }

  /** @return the total weight of all instances in the dataset */
  public double totalWeight() {
    return totalNegativeWeight() + totalPositiveWeight();
  }

  /** @return the total number of positive instances in the dataset (ignoring weight) */
  public long totalPositives() {
    return mPosCount;
  }

  /** @return the total number of negative instances in the dataset (ignoring weight) */
  public long totalNegatives() {
    return size() - totalPositives();
  }

  /** Equalize the total positive and negative weights of the dataset. */
  public void reweight() {
    // Includes a Laplace style correction to guard against potential of no positives or negatives
    final double targetWeightPerClass = 0.5 * (size() + 2);
    final double posWeight = targetWeightPerClass / (totalPositives() + 1);
    final double negWeight = targetWeightPerClass / (totalNegatives() + 1);
    mPosWeight = 0;
    mNegWeight = 0;
    for (final Instance i : getInstances()) {
      if (i.isPositive()) {
        i.setWeight(posWeight);
        mPosWeight += posWeight;
      } else {
        i.setWeight(negWeight);
        mNegWeight += negWeight;
      }
    }
  }

  /**
   * Randomly replace values within the dataset with the missing value
   * @param errorRate the probability of setting a value to missing
   */
  public void injectMissing(double errorRate) {
    injectMissing(errorRate, errorRate);
  }

  /**
   * Randomly replace values within the dataset with the missing value
   * @param posErrorRate the probability of setting a value to missing for positive instances
   * @param negErrorRate the probability of setting a value to missing for negative instances
   */
  public void injectMissing(double posErrorRate, double negErrorRate) {
    injectErrors(posErrorRate, negErrorRate, Double.NaN);
  }

  /**
   * Randomly replace values within the dataset with the specified value
   * @param posErrorRate the probability of setting a value to missing for positive instances
   * @param negErrorRate the probability of setting a value to missing for negative instances
   * @param newValue the value to inject
   */
  public void injectErrors(double posErrorRate, double negErrorRate, double newValue) {
    final PortableRandom rand = new PortableRandom(42);
    final int numAtts = getAttributes().length;
    for (Instance inst : getInstances()) {
      for (int i = 0; i < numAtts; ++i) {
        if (rand.nextDouble() < (inst.isPositive() ? posErrorRate : negErrorRate)) {
          inst.instance()[i] = newValue;
        }
      }
    }
  }

  /**
   * Get a count of the number of missing values for each attribute
   * @return the counts for each attribute
   */
  public long[] missingValueCounts() {
    final int numAtts = getAttributes().length;
    final long[] counts = new long[numAtts];
    for (Instance inst : getInstances()) {
      for (int attribute = 0; attribute < numAtts; attribute++) {
        if (Attribute.isMissingValue(inst.instance()[attribute])) {
          counts[attribute]++;
        }
      }
    }
    return counts;
  }

}
