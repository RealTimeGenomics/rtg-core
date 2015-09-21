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
package com.rtg.variant;

import java.io.IOException;
import java.util.Properties;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MathUtils;
import com.rtg.util.PropertiesUtils;
import com.rtg.util.Utils;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.machine.MachineType;

/**
 * A builder class for <code>MachineErrorParams</code>.
 */
@TestClass(value = { "com.rtg.variant.MachineErrorParamsTest" })
public class MachineErrorParamsBuilder {

  static final boolean CG_TRIM = true; //Boolean.valueOf(System.getProperty("cg-trim", "true"));

  static final double[] CG_DEFAULT_OVERLAP_DIST = MathUtils.renormalize(new int[]{0, 8, 84, 8, 0});
  static final double[] CG_DEFAULT_GAP_DIST = MathUtils.renormalize(new int[]{0, 27, 64, 9, 0});
  static final double[] CG_DEFAULT_SMALL_GAP_DIST = MathUtils.renormalize(new int[]{90, 7, 3, 0});
  static final double[] CG_DEFAULT_OVERLAP2_DIST = MathUtils.renormalize(new int[]{31, 58, 128, 376, 547, 86, 3, 0});

  /** The proportion of insertions that are homozygous */
  protected double mErrorMnpEventRate, mErrorInsEventRate, mErrorDelEventRate;

  protected double[] mErrorMnpDistribution, mErrorInsDistribution, mErrorDelDistribution;

  protected int[] mQualityCurve;

  // CG v1
  protected double[] mGapDistribution;
  protected double[] mSmallGapDistribution;
  protected double[] mOverlapDistribution;
  protected boolean mCGTrimOuterBases;

  // CG v2
  protected double[] mOverlapDistribution2;

  protected MachineType mMachine;


  /**
   * Creates a builder with initial default values from the
   * <code>human.properties</code> file.
   */
  public MachineErrorParamsBuilder() {
    try {
      errors("default");
    } catch (final Exception e) {
      throw new RuntimeException("Error reading default.properties", e);
    }
  }

  /**
   * Sets the name of the machine errors resource file to use. This reads the
   * resource file and sets all the probabilities.
   *
   * @param errors name of errors resource.
   * @throws InvalidParamsException if the resource file is invalid.
   * @throws IOException if the resource file cannot be read.
   */
  public MachineErrorParamsBuilder(String errors) throws IOException, InvalidParamsException {
    errors(errors);
  }

  /**
   * Construct builder from existing errors.
   *
   * @param errors name of errors resource.
   */
  public MachineErrorParamsBuilder(final AbstractMachineErrorParams errors) {
    mErrorMnpEventRate = errors.errorMnpEventRate();
    mErrorInsEventRate = errors.errorInsEventRate();
    mErrorDelEventRate = errors.errorDelEventRate();
    mErrorMnpDistribution = errors.errorMnpDistribution();
    mErrorInsDistribution = errors.errorInsDistribution();
    mErrorDelDistribution = errors.errorDelDistribution();
    mQualityCurve = errors.qualityCurve();
    mGapDistribution = errors.gapDistribution();
    mSmallGapDistribution = errors.smallGapDistribution();
    mOverlapDistribution = errors.overlapDistribution();
    mOverlapDistribution2 = errors.overlapDistribution2();
    mMachine = errors.machineType();
    mCGTrimOuterBases = errors.cgTrimOuterBases() && CG_TRIM;
  }

  /**
   * Sets the name of the machine errors resource file to use. This reads the
   * resource file and sets all the probabilities.
   *
   * @param errors name of errors resource. Default is default.
   * @return this builder, so calls can be chained.
   * @throws InvalidParamsException if the resource file is invalid.
   * @throws IOException if the resource file cannot be read.
   */
  public MachineErrorParamsBuilder errors(final String errors) throws InvalidParamsException, IOException {
    Diagnostic.developerLog("Loading machine errors for: " + errors);
    //System.out.println("Loading machine errors for: " + errors);
    final Properties pr = PropertiesUtils.getPriorsResource(errors, PropertiesUtils.PropertyType.ERROR_PROPERTY);
    mErrorMnpDistribution = GenomePriorParamsBuilder.parseDistribution(errors, pr, "error_mnp_distribution", 1);
    mErrorInsDistribution = GenomePriorParamsBuilder.parseDistribution(errors, pr, "error_ins_distribution", 1);
    mErrorDelDistribution = GenomePriorParamsBuilder.parseDistribution(errors, pr, "error_del_distribution", 1);
    mErrorMnpEventRate = GenomePriorParamsBuilder.getDouble(errors, pr, "error_mnp_event_rate");
    mErrorInsEventRate = GenomePriorParamsBuilder.getDouble(errors, pr, "error_ins_event_rate");
    mErrorDelEventRate = GenomePriorParamsBuilder.getDouble(errors, pr, "error_del_event_rate");

    // parse the optional quality calibration curve
    if (pr.containsKey(MachineErrorParams.QCALIB_KEY)) {
      final String curve = pr.getProperty(MachineErrorParams.QCALIB_KEY);
      final String[] nums = curve.split(", *");
      if (nums.length != MachineErrorParams.QUALITIES) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, curve, MachineErrorParams.QCALIB_KEY, errors);
      }
      mQualityCurve = new int[nums.length];
      for (int q = 0; q < nums.length; q++) {
        mQualityCurve[q] = Integer.parseInt(nums[q]);
      }
      checkQualityCurve(mQualityCurve); // check the values
    }
    if (pr.containsKey("machine_type")) {
      try {
        mMachine = MachineType.valueOf(pr.getProperty("machine_type"));
      } catch (final IllegalArgumentException e) {
        mMachine = null;
      }
    }
    mCGTrimOuterBases = tryGetBoolean(pr, "cg_trim_outer_bases", false) && CG_TRIM;
    mOverlapDistribution = getDistribution(pr, "overlap", CG_DEFAULT_OVERLAP_DIST, errors);
    mSmallGapDistribution = getDistribution(pr, "smallgap", CG_DEFAULT_SMALL_GAP_DIST, errors);
    mGapDistribution = getDistribution(pr, "gap", CG_DEFAULT_GAP_DIST, errors);
    mOverlapDistribution2 = getDistribution(pr, "overlap_2", CG_DEFAULT_OVERLAP2_DIST, errors);

    return this;
  }

  private double[] getDistribution(Properties pr, String key, double[] defaultDist, String errors) {
    final double[] dist;
    if (pr.containsKey(key)) {
      final String overlap = pr.getProperty(key);
      try {
        dist = MathUtils.renormalize(ArrayUtils.parseIntArray(overlap));
        if (dist.length < defaultDist.length) {
          throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, overlap, "overlap", errors);
        }
      } catch (final NumberFormatException e) {
        throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, overlap, "overlap", errors);
      }
    } else {
      dist = defaultDist;
    }
    return dist;
  }

  private static boolean tryGetBoolean(final Properties pr, final String key, final boolean dfault) {
    if (pr.containsKey(key)) {
      return Boolean.valueOf(pr.getProperty(key));
    } else {
      return dfault;
    }
  }

  static double parseDouble(final String prior, final String val, final String key) {
    final double ret;
    try {
      ret = Double.valueOf(val);
    } catch (final NumberFormatException e) {
      throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, val, key, prior);
    }
    if (ret < 0.0 || ret > 1.0) {
      throw new InvalidParamsException(ErrorType.PRIOR_KEY_VALUE_INVALID, val, key, prior);
    }
    return ret;
  }

  /**
   * Sets MNP event rate error prior.
   *
   * @param prior value.
   * @return this builder, so calls can be chained.
   */
  public MachineErrorParamsBuilder errorMnpEventRate(final double prior) {
    mErrorMnpEventRate = prior;
    return this;
  }

  /**
   * Sets insert error prior.
   *
   * @param prior value.
   * @return this builder, so calls can be chained.
   */
  public MachineErrorParamsBuilder errorInsEventRate(final double prior) {
    mErrorInsEventRate = prior;
    return this;
  }

  /**
   * Sets deletion error event prior.
   *
   * @param prior value.
   * @return this builder, so calls can be chained.
   */
  public MachineErrorParamsBuilder errorDelEventRate(final double prior) {
    mErrorDelEventRate = prior;
    return this;
  }

  /**
   * Set the length distribution of machine MNP errors. Note that
   * <code>lengths[0]</code> refers to inserts of length 1.
   *
   * @param lengths the probability of each length (1..). Must sum to 1.0.
   * @return this builder, so calls can be chained.
   * @throws InvalidParamsException if distribution is invalid.
   */
  public MachineErrorParamsBuilder errorMnpDistribution(final double[] lengths) {
    checkDistribution(lengths);
    mErrorMnpDistribution = new double[lengths.length + 1];
    System.arraycopy(lengths, 0, mErrorMnpDistribution, 1, lengths.length);
    return this;
  }

  /**
   * Set the length distribution of machine insertion errors. Note that
   * <code>lengths[0]</code> refers to inserts of length 1.
   *
   * @param lengths the probability of each length (1..). Must sum to 1.0.
   * @return this builder, so calls can be chained.
   * @throws InvalidParamsException if distribution is invalid.
   */
  public MachineErrorParamsBuilder errorInsDistribution(final double[] lengths) {
    checkDistribution(lengths);
    mErrorInsDistribution = new double[lengths.length + 1];
    System.arraycopy(lengths, 0, mErrorInsDistribution, 1, lengths.length);
    return this;
  }

  /**
   * Set the length distribution of machine deletion errors. Note that
   * <code>lengths[0]</code> refers to deletes of length 1.
   *
   * @param lengths the probability of each length (1..). Must sum to 1.0.
   * @return this builder, so calls can be chained.
   * @throws InvalidParamsException if distribution is invalid.
   */
  public MachineErrorParamsBuilder errorDelDistribution(final double[] lengths) {
    checkDistribution(lengths);
    mErrorDelDistribution = new double[lengths.length + 1];
    System.arraycopy(lengths, 0, mErrorDelDistribution, 1, lengths.length);
    return this;
  }

  /**
   * Set the flag that says whether to use CG Outer Base Trimming.
   *
   * @param cgTrimOuterBases flag - true if Trimming to be used
   * @return this builder, so calls can be chained.
   */
  public MachineErrorParamsBuilder cgTrimOuterBases(final boolean cgTrimOuterBases) {
    mCGTrimOuterBases = cgTrimOuterBases;
    return this;
  }

  /**
   * Set the machine type.
   *
   * @param machine the machine type.
   * @return this builder, so calls can be chained.
   */
  public MachineErrorParamsBuilder machine(final MachineType machine) {
    mMachine = machine;
    return this;
  }

  /**
   * Creates a MachineErrorParams using the current builder.
   *
   * @return the new parameters object.
   * @throws InvalidParamsException if any parameters are invalid.
   */
  public MachineErrorParams create() {
    checkProbability(mErrorInsEventRate);
    checkProbability(mErrorDelEventRate);
    checkDistribution(mErrorMnpDistribution);
    checkDistribution(mErrorInsDistribution);
    checkDistribution(mErrorDelDistribution);
    if (mQualityCurve != null) {
      checkQualityCurve(mQualityCurve);
    }
    checkDistribution(mGapDistribution);
    checkDistribution(mOverlapDistribution);
    checkDistribution(mSmallGapDistribution);
    checkDistribution(mOverlapDistribution2);
    return new MachineErrorParams(this);
  }

  private static void checkProbability(final double value) {
    if (value < 0.0 || value > 1.0) {
      throw new IllegalArgumentException("rate must be 0.0 .. 1.0, not " + Utils.realFormat(value, 6));
    }
  }

  static void checkDistribution(final double[] distrib) {
    double sum = 0.0;
    for (final double element : distrib) {
      checkProbability(element);
      sum += element;
    }
    if (Math.abs(sum - 1.0) > 0.0001) {
      throw new IllegalArgumentException("distribution must sum to 1.0, not " + sum);
    }
  }

  private static void checkQualityCurve(final int[] qualities) {
    if (qualities.length != MachineErrorParams.QUALITIES) {
      throw new InvalidParamsException(ErrorType.INVALID_QUALITY_LENGTH, MachineErrorParams.QCALIB_KEY);
    }
    for (int quality : qualities) {
      if (quality < 0 || quality >= MachineErrorParams.QUALITIES) {
        throw new InvalidParamsException(ErrorType.INVALID_QUALITY);
      }
    }
  }


}
