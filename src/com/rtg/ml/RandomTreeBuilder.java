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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Properties;

import com.rtg.launcher.GlobalFlags;
import com.rtg.ml.ZeroRBuilder.ZeroRClassifier;
import com.rtg.util.DoubleMultiSet;
import com.rtg.util.PortableRandom;
import com.rtg.util.Seedable;
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
  /** Property to set the random number seed */
  public static final String PROP_SEED = PROP_PREFIX + "seed";

  // Properties used during building, affect the resulting classifier
  /** Maximum depth of tree. 0 means no maximum */
  private int mMaxDepth = 0;
  /** Minimum number of instances at a leaf */
  private int mMinInstances = 10;
  /** Number of attributes used at each node. 0 means log base 2 (num attributes) + 1 */
  private int mNumAttributes = 0;
  /** Random number seed */
  private int mSeed = 42;


  private PredictClassifier mClassifier = null;
  private int mActualNumAttributes = 0;

  @Override
  public void build(Dataset dataset) {
    final PortableRandom random = new PortableRandom(mSeed);
    mActualNumAttributes = (mNumAttributes == 0) ? (int) (Math.log(dataset.getAttributes().length) / Math.log(2) + 1) : mNumAttributes;
    mClassifier = buildSubtree(random, dataset, 0);
  }

  private static final int POS = 0;
  private static final int NEG = 1;
  private static final int IN = 0;
  private static final int OUT = 1;

  private PredictClassifier buildSubtree(PortableRandom random, Dataset dataset, int currentDepth) {

    if ((dataset.size() < mMinInstances)
        || (mMaxDepth > 0 && currentDepth >= mMaxDepth)
        || (dataset.totalPositiveWeight() == 0)
        || (dataset.totalNegativeWeight() == 0)) {
      return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
    }

    // Choose set of attribute that the tree is allowed to select from
    BinarySplitter bestDirector = null;
    final double[][] priorDist = new double[2][2];
    priorDist[OUT][POS] = dataset.totalPositiveWeight();
    priorDist[OUT][NEG] = dataset.totalNegativeWeight();
    double bestEntropy = entropy(priorDist);
    //final long seed = random.getSeed();
    for (int attribute : getAttributes(random, dataset.getAttributes().length, mActualNumAttributes)) {
      // Evaluate each attribute for best split point
      final MlDataType dataType = dataset.getAttributes()[attribute].getDataType();
      if (dataType.isNumeric()) {

        // Create set of instances with non-missing attribute value, and initial counts
        final double[][] dist = new double[2][2];
        final ArrayList<Instance> sorted = new ArrayList<>();
        for (Instance inst : dataset.getInstances()) {
          if (!Attribute.isMissingValue(inst.instance()[attribute])) {
            sorted.add(inst);
            if (inst.isPositive()) {
              dist[OUT][POS] += inst.weight();
            } else {
              dist[OUT][NEG] += inst.weight();
            }
          }
        }
        final AttributeComparator comp = new AttributeComparator(attribute);
        Collections.sort(sorted, comp);

        // Scan through and find the best numeric split point
        double prevValue = Double.NaN;
        for (Instance inst : sorted) {
          final double currentValue = inst.instance()[attribute];
          if (!Attribute.isMissingValue(prevValue) && prevValue != currentValue) {
            // Evaluate gain
            final double entropy = entropy(dist);
            if (entropy < bestEntropy) {
              final double splitPoint = getSplitPoint(dataType, prevValue, currentValue);
              bestDirector = new BinarySplitter(dataset.getAttributes()[attribute].getName(), attribute, splitPoint, dataType);
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
//        final HashSet<Double> values = new HashSet<>();
//        values.addAll(posCounts.keySet());
//        values.addAll(negCounts.keySet());
        final double[][] dist = new double[2][2];
        // TODO is this the right thing to do with missing values
        final double posNonMissing = dataset.totalPositiveWeight() - posCounts.get(null);
        final double negNonMissing = dataset.totalNegativeWeight() - negCounts.get(null);
        final int nonimalSize = dataset.getAttributes()[attribute].nominalSize();
        for (int intKey = 0; intKey < nonimalSize; intKey++) {
          final Double key = (double) intKey;
          final double numpos = posCounts.get(key);
          final double numneg = negCounts.get(key);
          dist[IN][POS] = numpos;
          dist[IN][NEG] = numneg;
          dist[OUT][POS] = posNonMissing - numpos;
          dist[OUT][NEG] = negNonMissing - numneg;

          final double entropy = entropy(dist);
          if (entropy < bestEntropy) {
            bestDirector = new BinarySplitter(dataset.getAttributes()[attribute].getName(), attribute, key, dataType);
            bestEntropy = entropy;
          }
        }
      }
    }

    // Recurse for the chosen split point, or stop if no information gain
    if (bestDirector == null) {
      //System.out.println(dataset.totalPositives() + " " + dataset.totalNegatives() + " " + dataset.total());
      return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
    } else {

      final Dataset leftData = new Dataset(dataset.getAttributes());

      final double minWeight = filterInstances(bestDirector, dataset, leftData, null);

      // If all instances ever get filtered into the same branch then something has gone
      // wrong with the split point selection.  Ideally this should not happen, but
      // perhaps could with some combination of missing values.  In this situation
      // rather than further splitting, we just return a 0R on the input data.
      if (minWeight <= ZeroRClassifier.MINIMUM_WEIGHT) {
        Diagnostic.userLog("Unexpected empty branch during tree construction, using 0R instead of branching");
        return new ZeroRBuilder.ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
      }

      final double leftSize = leftData.size();
      final PredictClassifier left = buildSubtree(random, leftData, currentDepth + 1);
      // leftData can now be gc
      final Dataset rightData = new Dataset(dataset.getAttributes());

      filterInstances(bestDirector, dataset, null, rightData);
      final PredictClassifier right = buildSubtree(random, rightData, currentDepth + 1);
      final double leftFraction = leftSize / (leftSize + rightData.size());
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
   * @param left ends up with left instances (can be null if you don't want left)
   * @param right ends up with right instances (can be null if you don't want right)
   * @return weight of the smaller set (irrespective of whether it was null as a parameter)
   */
  private double filterInstances(BinarySplitter director, Dataset instances, Dataset left, Dataset right) {
    double leftWeight = 0;
    double rightWeight = 0;
    for (Instance inst : instances.getInstances()) {
      switch (director.split(inst.instance())) {
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
        if (GlobalFlags.getBooleanValue(GlobalFlags.AVR_TRAIN_ON_MISSING_VALUES)) {
          // Send instances with missing values down both branches with half weight
          final double halfWeight = 0.5 * inst.weight();
          leftWeight += halfWeight;
          rightWeight += halfWeight;
          if (left != null) {
            left.addInstance(inst.reweight(halfWeight));
          }
          if (right != null) {
            right.addInstance(inst.reweight(halfWeight));
          }
        }
        break;
      }
    }
    return Math.min(leftWeight, rightWeight);
  }

  static double entropy(double[][] dist) {
    return entropy(dist[IN][POS], dist[IN][NEG]) + entropy(dist[OUT][POS], dist[OUT][NEG]);
  }

  static double entropy(double pos, double neg) {
    final double tot = pos + neg;
    if (tot <= 0) {
      return 0;
    }
    return -(lnFunc(pos) + lnFunc(neg) - lnFunc(tot));
  }

  static double lnFunc(double p) {
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
    for (int i = 0; i < available.length; i++) {
      available[i] = i;
    }
    final int[] result = new int[subsetSize];
    int max = total;
    for (int i = 0; i < result.length; i++) {
      final int chosen = seed.nextInt(max);
      result[i] = available[chosen];
      max--;
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
    mMinInstances = Integer.parseInt(props.getProperty(PROP_MIN_INSTANCES, String.valueOf(10)));
    if (mMinInstances <= 0) {
      throw new IllegalArgumentException("Minimum number of instances must be greater than 1");
    }
    mMaxDepth = Integer.parseInt(props.getProperty(PROP_MAX_DEPTH, String.valueOf(10)));
    if (mMaxDepth < 0) {
      throw new IllegalArgumentException("Maximum tree depth cannot be negative");
    }
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
      @SuppressWarnings(value = {"unchecked", "rawtypes"})
      final int result = ((Comparable) a1).compareTo(a2);
      return result;
    }
  }
}
