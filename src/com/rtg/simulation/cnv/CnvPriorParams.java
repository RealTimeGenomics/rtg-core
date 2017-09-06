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
package com.rtg.simulation.cnv;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.Arrays;
import java.util.Properties;

import com.rtg.simulation.SimulationUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.ObjectParams;
import com.rtg.util.PropertiesUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.ErrorType;

/**
 */
public class CnvPriorParams extends ObjectParams {


  /**
   * Creates a CnvPriorParams builder.
   * @return the builder.
   */
  public static CnvPriorParamsBuilder builder() {
    return new CnvPriorParamsBuilder();
  }

  /**
   * A builder class for <code>CnvPriorParams</code>.
   */
  public static final class CnvPriorParamsBuilder {

    protected int[] mCopyRangeEnds = {1, 2, 5, 10, 20, 100};

    protected double[][] mCopyNumberDistribution;

    protected double[] mPowerLengthDistribution;

    protected double mProbDeletedOnOneStrand;

    protected double mProbDeletedOnBothStrands;
    private static final String[] MAGNITUTE_STR = {"1", "2", "3", "4", "5", "6", "7", "8"};

    /**
     * Creates a builder with initial default values from the <code>cnv-default.properties</code> file.
     */
    public CnvPriorParamsBuilder() {
      try {
        cnvpriors("cnv-default");
      } catch (final Exception e) {
        throw new RuntimeException("Error reading cnv-default.properties", e);
      }
    }

    /**
     * Sets the name of the priors resource to use.
     * This reads the resource file and sets all the probabilities.
     *
     * @param priorName name of priors resource. Default is human.
     * @return this builder, so calls can be chained.
     * @throws InvalidParamsException if the resource file is invalid.
     * @throws IOException if the resource file cannot be read.
     */
    public CnvPriorParamsBuilder cnvpriors(final String priorName) throws IOException {
      //TODO error messages change to cnv-priors
      final Properties pr = PropertiesUtils.getPriorsResource(priorName, PropertiesUtils.PropertyType.CNV_PROPERTY);

      mPowerLengthDistribution = parseDistribution(priorName, pr, "power_length_distribution", 0);

      mCopyRangeEnds = getIntegerArray(priorName, pr, "copy_range_ends");

      mCopyNumberDistribution = new double[MAGNITUTE_STR.length][];
      for (int i = 0; i < MAGNITUTE_STR.length; ++i) {
        mCopyNumberDistribution[i] = parseDistribution(priorName, pr,
            "copy_number_distribution_" + MAGNITUTE_STR[i], 0);
      }

      mProbDeletedOnOneStrand = getDouble(priorName, pr, "prob_deleted_on_one_strand");

      mProbDeletedOnBothStrands = getDouble(priorName, pr, "prob_deleted_on_both_strands");
      return this;
    }

    private double[] parseDistribution(final String prior, final Properties pr, final String key, final int start) {
      final String distrib = pr.getProperty(key);
      if (distrib == null) {
        throw new InvalidParamsException(ErrorType.PROPS_KEY_NOT_FOUND, key, prior);
      }
      final String[] words = distrib.split(",");
      final double[] result = new double[words.length + start];
      try {
        for (int i = 0; i < words.length; ++i) {
          result[i + start] = parseDouble(prior, words[i], key);
        }
        checkDistribution(result);
      } catch (final IllegalArgumentException e) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, distrib, key, prior);
      }
      return result;
    }

    private int[] getIntegerArray(final String prior, final Properties pr, final String key) {
      final String distrib = pr.getProperty(key);
      if (distrib == null) {
        throw new InvalidParamsException(ErrorType.PROPS_KEY_NOT_FOUND, key, prior);
      }
      final String[] words = distrib.split(",");
      final int[] result = new int[words.length];
      try {
        for (int i = 0; i < words.length; ++i) {
          result[i] = parseInteger(prior, words[i], key);
        }
      } catch (final IllegalArgumentException e) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, distrib, key, prior);
      }
      return result;
    }

    private static double getDouble(final String prior, final Properties pr, final String key) {
      final String val = pr.getProperty(key);
      //System.err.println("key=" + key + " val=" + val);
      if (val == null) {
        throw new InvalidParamsException(ErrorType.PROPS_KEY_NOT_FOUND, key, prior);
      }
      return parseDouble(prior, val, key);
    }

    private static double parseDouble(final String prior, final String val, final String key) {
      final double ret;
      try {
        ret = Double.parseDouble(val);
      } catch (final NumberFormatException e) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, val, key, prior);
      }
      if (ret < 0.0 || ret >= 1.0) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, val, key, prior);
      }
      return ret;
    }

    private static int parseInteger(final String prior, final String val, final String key) {
      final int ret;
      try {
        ret = Integer.parseInt(val);
      } catch (final NumberFormatException e) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, val, key, prior);
      }
      return ret;
    }

    /**
     * Sets copy number range distribution of large regions.
     *
     * @param rangeEnds list of range ends in order increasing (or equal) from left to right
     * @return this builder, so calls can be chained.
     * @throws InvalidParamsException if the distribution is invalid
     */
    public CnvPriorParamsBuilder copyRangeEnds(final int[] rangeEnds) {
      checkEndsList(rangeEnds);
      mCopyRangeEnds = Arrays.copyOf(rangeEnds, rangeEnds.length);
      return this;
    }

    /**
     * Creates a CnvPriorParams using the current builder.
     *
     * @return the new CnvPriorParams
     * @throws InvalidParamsException if any parameters are invalid.
     */
    public CnvPriorParams create() {
      // TODO: checks
      return new CnvPriorParams(this);
    }
  }

  protected static void checkDistribution(final double[] distrib) {
    double sum = 0.0;
    for (final double element : distrib) {
      checkProbability(element);
      sum += element;
    }
    if (Math.abs(sum - 1.0) > 0.0001) {
      throw new IllegalArgumentException("distribution must sum to 1.0, not " + sum);
    }
  }

  protected static void checkProbability(final double value) {
    if (value < 0.0 || value > 1.0
    ) {
      throw new IllegalArgumentException("rate must be 0.0 .. 1.0, not " + Utils.realFormat(value, 1));
    }
  }

  protected static void checkEndsList(final int[] endsList) {
    int last = 1;
    for (final int element : endsList) {
      if (element < last) {
        throw new IllegalArgumentException("ends list must contain values > 1 with from left to right increasing or equal values ");
      }
      last = element;
    }
  }

  private final int[] mCopyRangeEnds;
  private final double[][] mCopyNumberDistribution;
  private final double[] mPowerLengthDistribution;
  private final double[] mPowerLengthThresholds;
  private final double mProbDeletedOnBothStrands;
  private final double mProbDeletedOnOneStrand;
  private final double[][] mCopyNumberThresholds;

  /**
   * @param builder the builder object.
   */
  public CnvPriorParams(final CnvPriorParamsBuilder builder) {
    super();
    mCopyRangeEnds = builder.mCopyRangeEnds;
    mCopyNumberDistribution = builder.mCopyNumberDistribution;
    mPowerLengthDistribution = builder.mPowerLengthDistribution;
    mProbDeletedOnBothStrands = builder.mProbDeletedOnBothStrands;
    mProbDeletedOnOneStrand = builder.mProbDeletedOnOneStrand;

    mCopyNumberThresholds = new double[mCopyNumberDistribution.length][];
    for (int i = 0; i < mCopyNumberThresholds.length; ++i) {
      mCopyNumberThresholds[i] = SimulationUtils.cumulativeDistribution(mCopyNumberDistribution[i]);
    }
    mPowerLengthThresholds = SimulationUtils.cumulativeDistribution(mPowerLengthDistribution);
  }

  /**
   * Get the list of ends for the copy number ranges.
   * @return an integer array with the copy end values.
   */
  public int[] copyRangeEnds() {
    return mCopyRangeEnds;
  }

  /**
   * Get  distribution for large regions.
   *
   * @return a (read-only) array of doubles that sums to 1.0.
   */
  public double[][] copyNumberDistribution() {
    return mCopyNumberDistribution;
  }

  /**
   * Get thresholds of  distribution for small regions.
   *
   * @return a (read-only) array of doubles in ascending order up to 1.0.
   */
  public double[][] copyNumberThresholds() {
    return mCopyNumberThresholds;
  }

  /**
   * Get distribution of the powers of the region lengths.
   *
   * @return a (read-only) array of doubles that sums to 1.0.
   */
  public double[] powerLengthDistribution() {
    return mPowerLengthDistribution;
  }

  /**
   * Get thresholds of the  of the powers of the region lengths.
   *
   * @return a (read-only) array of doubles in ascending order up to 1.0.
   */
  public double[] powerLengthThresholds() {
    return mPowerLengthThresholds;
  }

  /**
   * Gets the probability for the region deleted on both strands
   *
   * @return the probability for the region deleted on both strands
   */
  public double probDeletedOnBothStrands() {
    return mProbDeletedOnBothStrands;
  }

  /**
   * Gets the probability for the region deleted on one strand
   *
   * @return the probability for the region deleted on one strand
   */
  public double probDeletedOnOneStrand() {
    return mProbDeletedOnOneStrand;
  }

  @Override
  public String toString() {
    return "    "
    + "probability one delete = " + probDeletedOnOneStrand() + LS
    + "probability both delete = " + probDeletedOnBothStrands() + LS;
  }

}
