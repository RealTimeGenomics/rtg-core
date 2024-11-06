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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import java.util.Properties;

import com.rtg.ml.ZeroRBuilder.ZeroRClassifier;
import com.rtg.util.DoubleMultiSet;
import com.rtg.util.PortableRandom;
import com.rtg.util.Seedable;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Build a decision tree with the set of attributes available at each node being a random subset of all attributes.
 *
 */
public class RandomTreeBuilder implements BuildClassifier, Seedable {

  /** All our properties start with this */
  public static final String PROP_PREFIX = "randomtree.";

  /** Property to set the maximum depth of tree that can be constructed */
  public static final String PROP_MAX_DEPTH = PROP_PREFIX + "max-depth";
  /** Property to set the minimum number of instances before allowing split */
  public static final String PROP_MIN_INSTANCES = PROP_PREFIX + "min-instances";
  /** Property to set the number of attributes to consider at each node */
  public static final String PROP_NUM_ATTRIBUTES = PROP_PREFIX + "num-attributes";
  /** Property to set additional split cost */
  public static final String PROP_SPLIT_COST = PROP_PREFIX + "split-cost";
  /** Property to set how instances with missing values should be sent into subtrees */
  public static final String PROP_PROPAGATE_MISSING = PROP_PREFIX + "propagate-missing";
  /** Property to set whether to use new entropy calculation which includes missing values */
  public static final String PROP_ENTROPY_MISSING = PROP_PREFIX + "entropy-missing";
  /** Property to set whether to consider a split point that uses a missing vs non-missing distinction */
  public static final String PROP_SPLIT_MISSING = PROP_PREFIX + "split-missing";
  /** Property to set the random number seed */
  public static final String PROP_SEED = PROP_PREFIX + "seed";

  private enum PropagateMissingType {
    /**Don't send missing value instances into sub-trees */
    FALSE,
    /** Send missing value instances into both sub-trees, weighted according to fraction */
    BOTH,
    /** Send missing value instances into a randomly selected sub-tree, according to fraction */
    RANDOM,
  }

  // Properties used during building, affect the resulting classifier
  /** Maximum depth of tree. 0 means no maximum */
  private int mMaxDepth = 0;
  /** Minimum number of instances at a leaf */
  private int mMinInstances = 10;
  /** Number of attributes used at each node. 0 means log base 2 (num attributes) + 1 */
  private int mNumAttributes = 0;
  /** Optional additional cost to splitting */
  private double mSplitCost = 0;
  /** Random number seed */
  private int mSeed = 42;
  /** Sub-tree behaviour for instances with missing values */
  private PropagateMissingType mPropagateMissing = PropagateMissingType.FALSE;
  /** Use entropy calculation which includes missing values */
  private boolean mEntropyMissing = false;
  /** Consider a missing vs non-missing split */
  private boolean mSplitMissing = false;

  
  private PredictClassifier mClassifier = null;
  private int mActualNumAttributes = 0;

  @Override
  public void build(Dataset dataset) {
    final PortableRandom random = new PortableRandom(mSeed);
    mActualNumAttributes = (mNumAttributes == 0) ? (int) (Math.log(dataset.getAttributes().length) / Math.log(2) + 1) : mNumAttributes;
    Diagnostic.userLog(toString());
    mClassifier = buildSubtree(random, dataset, 0);
  }

  private static final int POS = 0;
  private static final int NEG = 1;
  private static final int IN = 0;
  private static final int OUT = 1;

  private PredictClassifier buildSubtree(PortableRandom random, Dataset dataset, int currentDepth) {

    if ((dataset.totalWeight() < mMinInstances)
        || (mMaxDepth > 0 && currentDepth >= mMaxDepth)
        || (dataset.totalPositiveWeight() == 0)
        || (dataset.totalNegativeWeight() == 0)) {
      return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
    }

    // Choose set of attribute that the tree is allowed to select from
    BinarySplitter bestDirector = null;
    double bestFrac = Double.NaN;
    final double[][] priorDist = new double[2][2];
    priorDist[OUT][POS] = dataset.totalPositiveWeight();
    priorDist[OUT][NEG] = dataset.totalNegativeWeight();
    double bestEntropy = entropy(priorDist) - mSplitCost;
    //final long seed = random.getSeed();
    for (int attribute : getAttributes(random, dataset.getAttributes().length, mActualNumAttributes)) {
      // Evaluate each attribute for best split point
      final Attribute att = dataset.getAttributes()[attribute];
      final MlDataType dataType = att.getDataType();
      if (dataType.isNumeric()) {

        // Create set of instances with non-missing attribute value, and initial counts
        final double[][] dist = new double[2][2];
        final ArrayList<Instance> sorted = new ArrayList<>();
        for (Instance inst : dataset.getInstances()) {
          if (!Attribute.isMissingValue(inst.instance()[attribute])) {
            sorted.add(inst);
            dist[OUT][inst.isPositive() ? POS : NEG] += inst.weight();
          }
        }
        final double missingPos = dataset.totalPositiveWeight() - dist[OUT][POS];
        final double missingNeg = dataset.totalNegativeWeight() - dist[OUT][NEG];
        if (mSplitMissing) {
          assert mEntropyMissing;
          final double entropy = entropy(dist) + entropy(missingPos, missingNeg);
          if (entropy < bestEntropy) {
            bestDirector = new BinarySplitter(att.getName(), attribute, Double.NaN, dataType);
            bestEntropy = entropy;
          }
        }

        sorted.sort(new AttributeComparator(attribute));

        // Scan through and find the best numeric split point
        double prevValue = Double.NaN;
        for (Instance inst : sorted) {
          final double currentValue = inst.instance()[attribute];
          if (prevValue != currentValue && !Attribute.isMissingValue(prevValue)) {
            // Evaluate gain
            final double entropy = mEntropyMissing ? entropy(dist, missingPos, missingNeg) : entropy(dist);
            if (entropy < bestEntropy) {
              final double splitPoint = getSplitPoint(dataType, prevValue, currentValue);
              bestDirector = new BinarySplitter(att.getName(), attribute, splitPoint, dataType);
              bestFrac = in(dist) / total(dist);

              assert BinarySplitter.Direction.LEFT == bestDirector.split(prevValue);
              assert BinarySplitter.Direction.LEFT == bestDirector.split(splitPoint);
              assert BinarySplitter.Direction.RIGHT == bestDirector.split(currentValue);
              bestEntropy = entropy;
            }
          }
          // Update dist
          if (inst.isPositive()) {
            dist[IN][POS] += inst.weight();
            dist[OUT][POS] -= inst.weight();
          } else {
            dist[IN][NEG] += inst.weight();
            dist[OUT][NEG] -= inst.weight();
          }
          prevValue = currentValue;
        }

      } else { // Nominal attributes
        // Collect distribution for all distinct values
        final DoubleMultiSet<Double> posCounts = new DoubleMultiSet<>();
        final DoubleMultiSet<Double> negCounts = new DoubleMultiSet<>();
        for (Instance inst : dataset.getInstances()) {
          final double attValue = inst.instance()[attribute];
          if (inst.isPositive()) {
            posCounts.add(Attribute.isMissingValue(attValue) ? null : attValue, inst.weight());
          } else {
            negCounts.add(Attribute.isMissingValue(attValue) ? null : attValue, inst.weight());
          }
        }
        final double[][] dist = new double[2][2];
        final double missingPos = posCounts.get(null);
        final double missingNeg = negCounts.get(null);
        if (mSplitMissing) {
          assert mEntropyMissing;
          final double entropy = entropy(dist) + entropy(missingPos, missingNeg);
          if (entropy < bestEntropy) {
            bestDirector = new BinarySplitter(att.getName(), attribute, Double.NaN, dataType);
            bestEntropy = entropy;
          }
        }
        final double posNonMissing = dataset.totalPositiveWeight() - missingPos;
        final double negNonMissing = dataset.totalNegativeWeight() - missingNeg;
        final int nonimalSize = att.nominalSize();
        double bestKey = Double.NaN;
        for (int intKey = 0; intKey < nonimalSize; ++intKey) {
          final Double key = (double) intKey;
          final double numpos = posCounts.get(key);
          final double numneg = negCounts.get(key);
          dist[IN][POS] = numpos;
          dist[IN][NEG] = numneg;
          dist[OUT][POS] = posNonMissing - numpos;
          dist[OUT][NEG] = negNonMissing - numneg;

          final double entropy = mEntropyMissing ? entropy(dist, missingPos, missingNeg) : entropy(dist);
          if ((entropy < bestEntropy)
            || (entropy == bestEntropy && !Double.isNaN(bestKey) && att.compare(key, bestKey) < 0)) {
            bestDirector = new BinarySplitter(att.getName(), attribute, key, dataType);
            bestFrac = in(dist) / total(dist);
            bestEntropy = entropy;
            bestKey = key;
          }
        }
      }
    }

    // Recurse for the chosen split point, or stop if no information gain
    if (bestDirector == null) {
      return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());

    } else {
      final long seed = random.nextLong(); // We need to use this same seed for both times we run the filterInstances.
      Dataset subData = new Dataset(dataset.getAttributes());
      // bestFrac = 0.5; // Equivalent to old behaviour
      final double minWeight = filterInstances(bestDirector, bestFrac, dataset, subData, null, mPropagateMissing == PropagateMissingType.RANDOM ? new PortableRandom(seed) : null);
      // If all instances ever get filtered into the same branch then something has gone
      // wrong with the split point selection.  Ideally this should not happen, but
      // perhaps could with some combination of missing values.  In this situation
      // rather than further splitting, we just return a 0R on the input data.
      if (minWeight <= ZeroRClassifier.MINIMUM_WEIGHT) {
        Diagnostic.userLog("Unexpected empty branch during tree construction, using 0R instead of branching");
        return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
      }
      final double leftSize = subData.totalWeight();
      final PredictClassifier left = buildSubtree(random, subData, currentDepth + 1);

      subData = new Dataset(dataset.getAttributes());
      filterInstances(bestDirector, bestFrac, dataset, null, subData, mPropagateMissing == PropagateMissingType.RANDOM ? new PortableRandom(seed) : null);
      final PredictClassifier right = buildSubtree(random, subData, currentDepth + 1);

      final double leftFraction = leftSize / (leftSize + subData.totalWeight());
      return new BinaryTreeClassifier(bestDirector, left, right,  leftFraction);
    }
  }

  static double getSplitPoint(MlDataType dataType, double prevValue, double currentValue) {
    // WARNING: The treatment of integer attributes here is not type safe and
    // makes assumptions about how they are encoded as doubles.
    final double splitPoint;
    if (dataType == MlDataType.DOUBLE) {
      // Need to be careful with rounding here.  The "average" of two doubles can be the
      // larger of the two values, we need to ensure we always pick the smaller of the
      // two in this situation.
      final double s = (currentValue + prevValue) / 2;
      splitPoint = s == currentValue ? prevValue : s;
    } else if (dataType == MlDataType.INTEGER) {
      final int intSplit = ((int) currentValue + (int) prevValue) / 2;
      splitPoint = (double) intSplit;
    } else {
      splitPoint = prevValue;
    }
    return splitPoint;
  }

  /**
   * Split the instances into left and right according to the director.  Returns the weight
   * of the smaller of the two sides (useful for aborting computation).  Supports two pass
   * approach by setting left or right to null.
   * @param director how to choose
   * @param leftFraction fraction of missing value instances to assign to the left subtree
   * @param left ends up with left instances (can be null if you don't want left)
   * @param right ends up with right instances (can be null if you don't want right)
   * @param random if non-null, use stochastic selection of subtree for missing values
   * @return weight of the smaller set (irrespective of whether it was null as a parameter)
   */
  private double filterInstances(BinarySplitter director, double leftFraction, Dataset instances, Dataset left, Dataset right, PortableRandom random) {
    double leftWeight = 0;
    double rightWeight = 0;
    for (Instance inst : instances.getInstances()) {
      BinarySplitter.Direction d = director.split(inst.instance());
      if (d == BinarySplitter.Direction.MISSING && random != null) {
        if (Double.isNaN(leftFraction)) {
          throw new IllegalStateException("leftFraction should be set");
        }
        d = random.nextDouble() <= leftFraction ? BinarySplitter.Direction.LEFT : BinarySplitter.Direction.RIGHT;
      }
      switch (d) {
      case LEFT:
        leftWeight += inst.weight();
        if (left != null) {
          left.addInstance(inst);
        }
        break;
      case RIGHT:
        rightWeight += inst.weight();
        if (right != null) {
          right.addInstance(inst);
        }
        break;
      case MISSING:
      default:
        if (mPropagateMissing == PropagateMissingType.BOTH) {
          assert !Double.isNaN(leftFraction);
          // Send instances with missing values down both branches with proportional weight
          final double lWeight = leftFraction * inst.weight();
          final double rWeight = (1.0 - leftFraction) * inst.weight();
          leftWeight += lWeight;
          rightWeight += rWeight;
          if (left != null) {
            left.addInstance(inst.reweight(lWeight));
          }
          if (right != null) {
            right.addInstance(inst.reweight(rWeight));
          }
        }
        break;
      }
    }
    return Math.min(leftWeight, rightWeight);
  }

  private static double entropy(double[][] dist) {
    return entropy(dist, 0, 0);
  }

  private static double entropy(double[][] dist, double missingPos, double missingNeg) {
    final double fraction = in(dist) / total(dist);
    return entropy(dist[IN][POS] + fraction * missingPos, dist[IN][NEG] + fraction * missingNeg)
      + entropy(dist[OUT][POS] + (1 - fraction) * missingPos, dist[OUT][NEG] + (1 - fraction) * missingNeg);
  }

  private static double in(double[][] dist) {
    return ArrayUtils.sum(dist[IN]);
  }
  private static double out(double[][] dist) {
    return ArrayUtils.sum(dist[OUT]);
  }
  private static double total(double[][] dist) {
    return in(dist) + out(dist);
  }
  private static double entropy(double pos, double neg) {
    final double tot = pos + neg;
    if (tot <= 0) {
      return 0;
    }
    return -(lnFunc(pos) + lnFunc(neg) - lnFunc(tot));
  }

  private static double lnFunc(double p) {
    if (p <= 0) {
      return 0;
    }
    return p * Math.log(p);
  }

  /**
   * Choose a random set of attributes
   * @param seed random number seed
   * @param total total number of attributes
   * @param subsetSize number of attributes to select
   * @return an array of attribute ids
   */
  static int[] getAttributes(PortableRandom seed, int total, int subsetSize) {
    final int[] available = new int[total];
    for (int i = 0; i < available.length; ++i) {
      available[i] = i;
    }
    if (subsetSize >= total) {
      return available;
    }
    final int[] result = new int[subsetSize];
    int max = total;
    for (int i = 0; i < result.length; ++i) {
      final int chosen = seed.nextInt(max);
      result[i] = available[chosen];
      --max;
      available[chosen] = available[max];
    }
    return result;
  }

  @Override
  public PredictClassifier getClassifier() {
    return mClassifier;
  }

  @Override
  public void setProperties(Properties props) {
    setSeed(Integer.parseInt(props.getProperty(PROP_SEED, "42")));
    mNumAttributes = Integer.parseInt(props.getProperty(PROP_NUM_ATTRIBUTES, "0"));
    if (mNumAttributes < 0) {
      throw new IllegalArgumentException("Number of attributes can not be negative");
    }
    mMinInstances = Integer.parseInt(props.getProperty(PROP_MIN_INSTANCES, "10"));
    if (mMinInstances <= 0) {
      throw new IllegalArgumentException("Minimum number of instances must be greater than 1");
    }
    mMaxDepth = Integer.parseInt(props.getProperty(PROP_MAX_DEPTH, "10"));
    if (mMaxDepth < 0) {
      throw new IllegalArgumentException("Maximum tree depth cannot be negative");
    }
    mSplitCost = Double.parseDouble(props.getProperty(PROP_SPLIT_COST, "0"));
    mEntropyMissing = Boolean.parseBoolean(props.getProperty(PROP_ENTROPY_MISSING, Boolean.TRUE.toString()));
    mPropagateMissing = PropagateMissingType.valueOf(props.getProperty(PROP_PROPAGATE_MISSING, PropagateMissingType.RANDOM.toString()).toUpperCase(Locale.ROOT));
    mSplitMissing = Boolean.parseBoolean(props.getProperty(PROP_SPLIT_MISSING, Boolean.TRUE.toString()));
  }

  @Override
  public String toString() {
    return "RandomTreeBuilder:"
      + " " + PROP_SEED + "=" + mSeed
      + " " + PROP_NUM_ATTRIBUTES + "=" + mNumAttributes
      + " " + PROP_MIN_INSTANCES + "=" + mMinInstances
      + " " + PROP_MAX_DEPTH + "=" + mMaxDepth
      + " " + PROP_SPLIT_COST + "=" + mSplitCost
      + " " + PROP_ENTROPY_MISSING + "=" + mEntropyMissing
      + " " + PROP_PROPAGATE_MISSING + "=" + mPropagateMissing.name().toLowerCase(Locale.ROOT)
      + " " + PROP_SPLIT_MISSING + "=" + mSplitMissing;
  }

  @Override
  public int getSeed() {
    return mSeed;
  }

  @Override
  public void setSeed(int seed) {
    mSeed = seed;
  }

  static class AttributeComparator implements Comparator<Instance>, Serializable {
    private final int mAttributeIndex;

    AttributeComparator(int attributeIndex) {
      mAttributeIndex = attributeIndex;
    }

    @Override
    @SuppressWarnings(value = {"unchecked", "rawtypes"})
    public int compare(Instance o1, Instance o2) {
      final double a1 = o1.instance()[mAttributeIndex];
      final double a2 = o2.instance()[mAttributeIndex];
      if (a1 == a2) {
        return 0;
      }
      if (Attribute.isMissingValue(a1)) {
        return -1;
      } else if (Attribute.isMissingValue(a2)) {
        return 1;
      }
      return ((Comparable) a1).compareTo(a2);
    }
  }
}
