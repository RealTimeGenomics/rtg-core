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
package com.rtg.variant.util;


import java.util.Arrays;
import java.util.Locale;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.Arm;
import com.rtg.sam.MateInfo;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.MapInfo;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public final class VariantUtils {

  private static final String MAX_POSTERIOR = "1.79E308";
  //  private static final byte[] MAX_POSTERIOR_BYTES = MAX_POSTERIOR.getBytes();

  /** The value for the variant sample hypothesis name separator */
  public static final char COLON = ':';

  private VariantUtils() { } //prevent instantiation

  /**
   * The maximum value allowed in sam file quality fields
   */
  private static final int MAX_PHRED = 5000;
  private static final double[] PHRED_PROB;
  static final double BASE_LOG = Math.log1p(1.0);
  static final int BITS = 24;
  /** The minimum value of X that will return a useful value from <code>log1pExpApproximation</code> */
  public static final double MIN_X = 0.0001;
  static final int START_INDEX = doubleToIndex(MIN_X);
  static final double MAX_X;
  static final double[] VALUES;
  static final int SIZE;

  static {
    PHRED_PROB = new double[MAX_PHRED];
    for (int i = 0; i < MAX_PHRED; ++i) {
      PHRED_PROB[i] = Math.pow(10.0, -i / 10.0);
    }

    final int lastIndex;
    double x = MIN_X;
    for (; true; x += 0.01) {
      final double f = logExpFunction(x);
      if (f < MIN_X / 4.0) {
        MAX_X = x;
        lastIndex = doubleToIndex(x);
        SIZE = lastIndex - START_INDEX + 1;
        break;
      }
    }
    VALUES = new double[SIZE];
    for (int i = START_INDEX; i <= lastIndex; ++i) {
      final double u1 = indexToDouble(i);
      final double u2 = indexToDouble(i + 1);
      final double u = (u1 + u2) / 2.0;
      VALUES[i - START_INDEX] = logExpFunction(u);
    }
  }

  /**
   * Convert from a phred score to a probability.
   * @param phred phred formatted score to be converted.
   * @return probability of an error given the phred score.
   */
  public static double phredToProb(final int phred) {
    if (phred >= MAX_PHRED) {
      // Hopefully a rare case ...
      return Math.pow(10.0,  -phred / 10.0);
    }
    return PHRED_PROB[phred];
  }

  /**
   * Compute sum of two logs that is: log(exp(x) + exp(y)). Care is taken to
   * avoid over/under flow in the exp calculation. Can probably be made faster.
   *
   * @param x first log to be added.
   * @param y second log to be added.
   * @return <code>log(exp(x) + exp(y))</code>.
   */
  public static double logSum(final double x, final double y) {
    assert !Double.isNaN(x);
    assert !Double.isNaN(y);
    final double xp;
    final double yp;
    if (x < y) {
      xp = y;
      yp = x;
    } else {
      xp = x;
      yp = y;
    }
    if (yp == Double.NEGATIVE_INFINITY) {
      return xp;
    }
    if (xp == Double.POSITIVE_INFINITY) {
      return Double.POSITIVE_INFINITY;
    }
    final double z = yp - xp;
    final double res = xp + Math.log1p(Math.exp(z));
    assert !Double.isNaN(res) : xp + ":" + yp;
    return res;
  }

  /**
   * Compute subtraction of two logs that is: log(exp(x) - exp(y)). Care is taken to
   * avoid over/under flow in the exp calculation. Can probably be made faster.
   *
   * @param x first log value.
   * @param y second log to be subtracted.
   * @return <code>log(exp(x) - exp(y))</code>.
   */
  public static double logSubtract(final double x, final double y) {
    assert !Double.isNaN(x);
    assert !Double.isNaN(y);
    if (x < y) {
      return Double.NaN;
    }
    if (y == Double.NEGATIVE_INFINITY) {
      return x;
    }
    if (x == Double.POSITIVE_INFINITY) {
      return (y == Double.POSITIVE_INFINITY) ? Double.NaN : Double.POSITIVE_INFINITY;
    }
    final double z = y - x;
    final double res = x + Math.log1p(-Math.exp(z));
    assert !Double.isNaN(res) : x + ":" + y;
    return res;
  }


  /**
   * Format a posterior score in a common format.
   * @param q posterior score (in ln units) to be written.
   * @return the formatted posterior score
   */
  public static String formatPosterior(double q) {
    final double ql = q / MathUtils.LOG_10;
    //System.err.println("fpos: " + ql);
    if (Double.isInfinite(ql) && ql > 0.0) {
      return MAX_POSTERIOR;
    } else {
      return Utils.realFormat(ql, 1);
    }
  }

  /**
   * Compute sum of two logs that is: log(exp(x) + exp(y)). Care is taken to
   * avoid over/under flow in the exp calculation. Approximation only.
   *
   * @param x first log to be added.
   * @param y second log to be added.
   * @return <code>log(exp(x) + exp(y))</code>.
   */
  public static double logSumApproximation(final double x, final double y) {
    assert !Double.isNaN(x);
    assert !Double.isNaN(y);
    final double xp;
    final double yp;
    if (x < y) {
      xp = y;
      yp = x;
    } else {
      xp = x;
      yp = y;
    }
    if (yp == Double.NEGATIVE_INFINITY) {
      return xp;
    }
    if (xp == Double.POSITIVE_INFINITY) {
      return Double.POSITIVE_INFINITY;
    }
    final double z = xp - yp;
    final double res = xp + log1pExpApproximation(z);
    assert !Double.isNaN(res) : xp + ":" + yp;
    return res;
  }

  /**
   * Fast approximate calculation of <code>ln(1.0 + exp(-x))</code>
   * @param x argument for function.
   * @return result of function. <code>:p</code>
   */
  public static double log1pExpApproximation(final double x) {
    assert x >= 0.0;
    if (x < MIN_X) {
      return BASE_LOG;
    }
    if (x > MAX_X) {
      return 0;
    }
    final int i = doubleToIndex(x) - START_INDEX;
    if (i < 0 || i >= VALUES.length) {
      throw new RuntimeException("log1pExpApproximation i=" + i + " x=" + x + " START_INDEX=" + START_INDEX + " doubleToIndex=" + doubleToIndex(x) + " MIN_X=" + MIN_X + " MAX_X=" + MAX_X);
    }
    return VALUES[i];
  }

  static double indexToDouble(final int index) {
    return Double.longBitsToDouble(((long) index) << (64 - BITS));
  }

  static int doubleToIndex(final double x) {
    final long j = Double.doubleToRawLongBits(x);
    final long k = j >>> (64 - BITS);
    assert k <= Integer.MAX_VALUE;
    return (int) k;
  }

  /**
   * Compute <code>log(1+exp(-x))</code>.
   * @param x argument for function.
   * @return result of function.
   */
  public static double logExpFunction(final double x) {
    return Math.log1p(Math.exp(-x));
  }

  /**
   * A Main Method for Testing
   * @param args the arguments
   */
  public static void main(final String[] args) {
    System.out.println("BITS=" + BITS);
    System.out.println("MIN_X=" + MIN_X);
    System.out.println("MAX_X=" + MAX_X);
    System.out.println("START_INDEX=" + START_INDEX);
    System.out.println("SIZE=" + SIZE);
    double worstX = 0.0;
    double worstDelta = 0.0;
    double worstError = 0.0;
    for (double x = 0.0; x < 50.0; x += 0.001) {
      final double y = log1pExpApproximation(x);
      final double e = logExpFunction(x);
      final double delta = y - e;
      final double error = Math.abs(delta);
      if (error > worstError) {
        worstX = x;
        worstDelta = delta;
        worstError = error;
      }
      if (error > MIN_X) {
        System.out.println("x=" + x + " delta=" + delta);
      }
    }
    System.out.println("worstDelta=" + worstDelta + " at x=" + worstX);
  }

  /**
   * @param rec the SAM record to be read from.
   * @param params the parameters
   * @param <T> a MateInfo/MapInfo implementer
   * @return the score (Phred format)
   */
  public static <T extends MateInfo & MapInfo> int readScoreFromAlignmentRecord(final T rec, final VariantParams params) {
    final int readScore;
    final boolean isMated = rec.isMated();
    final int rs = rec.getMappingQuality();
    if (rs == 255 || params.ignoreReadQualities()) {
      final int nh = rec.getNHOrIH();
      if (nh == -1 || nh == 1) {
        readScore = isMated ? params.matedReadDefault() : params.unmatedReadDefault();
      } else if (nh == 0) {
        return 0;
      } else {
        readScore = phredFromN(nh);
      }
    } else {
      readScore = rs;
    }
    final int max = isMated ? params.matedReadMax() : params.unmatedReadMax();
    return Math.min(max, readScore);
  }

  /**
   * Dump string form of machine errors
   * @param p the machine errors
   * @return string of machine errors
   */
  public static String dumpMachineErrors(AbstractMachineErrorParams p) {
    final StringBuilder sb = new StringBuilder();
    sb.append("errorSnpRate: ").append(p.errorSnpRate()).append(StringUtils.LS);
    sb.append("errorMnpEventRate: ").append(p.errorMnpEventRate()).append(StringUtils.LS);
    sb.append("errorMnpDistribution: ").append(Arrays.toString(p.errorMnpDistribution())).append(StringUtils.LS);
    sb.append("errorDelBaseRate: ").append(p.errorDelBaseRate()).append(StringUtils.LS);
    sb.append("errorDelEventRate: ").append(p.errorDelEventRate()).append(StringUtils.LS);
    sb.append("errorDelDistribution: ").append(Arrays.toString(p.errorDelDistribution())).append(StringUtils.LS);
    sb.append("errorInsBaseRate: ").append(p.errorInsBaseRate()).append(StringUtils.LS);
    sb.append("errorInsEventRate: ").append(p.errorInsEventRate()).append(StringUtils.LS);
    sb.append("errorInsDistribution: ").append(Arrays.toString(p.errorInsDistribution())).append(StringUtils.LS);
    sb.append("machineType: ").append(p.machineType()).append(StringUtils.LS);
    sb.append("isCG: ").append(p.isCG()).append(StringUtils.LS);
    sb.append("overlapDistribution: ").append(Arrays.toString(p.overlapDistribution())).append(StringUtils.LS);
    sb.append("gapDistribution: ").append(Arrays.toString(p.gapDistribution())).append(StringUtils.LS);
    sb.append("smallGapDistribution: ").append(Arrays.toString(p.smallGapDistribution())).append(StringUtils.LS);
    sb.append("overlapDistribution2: ").append(Arrays.toString(p.overlapDistribution2())).append(StringUtils.LS);
    sb.append("quality curve: ");
    for (Arm arm : Arm.values()) {
      sb.append(arm).append("=");
      for (byte q = 0; q < 64; ++q) {
        sb.append(p.getScaledPhred(q, 0, arm)).append(",");
      }
      sb.append(StringUtils.LS);
    }
    sb.append(StringUtils.LS);
    return sb.toString();
  }

  private static String toString(double[] arr) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < arr.length; ++i) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(String.format(Locale.ROOT, "%.3g", arr[i]));
    }
    return sb.toString();
  }

  /**
   * Dump string form of machine errors
   * @param p the machine errors
   * @return string of machine errors
   */
  public static String toMachineErrorProperties(AbstractMachineErrorParams p) {
    final StringBuilder sb = new StringBuilder();
    if (p.machineType() != null) {
      sb.append("# machineType: ").append(p.machineType()).append(StringUtils.LS);
    }
    sb.append("# isCG: ").append(p.isCG()).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    sb.append(String.format(Locale.ROOT, "# derived: error_snp_rate = %.6f" + StringUtils.LS, p.errorSnpRate()));
    sb.append(String.format(Locale.ROOT, "error_mnp_event_rate = %.6f" + StringUtils.LS, p.errorMnpEventRate()));
    sb.append("error_mnp_distribution = ").append(toString(p.errorMnpDistribution())).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    sb.append(String.format(Locale.ROOT, "# derived: error_del_base_rate = %.6f" + StringUtils.LS, p.errorDelBaseRate()));
    sb.append(String.format(Locale.ROOT, "error_del_event_rate = %.6f" + StringUtils.LS, p.errorDelEventRate()));
    sb.append("error_del_distribution = ").append(toString(p.errorDelDistribution())).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    sb.append(String.format(Locale.ROOT, "# derived: error_ins_base_rate = %.6f" + StringUtils.LS, p.errorInsBaseRate()));
    sb.append(String.format(Locale.ROOT, "error_ins_event_rate = %.6f" + StringUtils.LS, p.errorInsEventRate()));
    sb.append("error_ins_distribution = ").append(toString(p.errorInsDistribution())).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    sb.append("quality_curve = ");
    for (Arm arm : Arm.values()) {
      sb.append(arm).append("=");
      for (byte q = 0; q < 64; ++q) {
        if (q != 0) {
          sb.append(", ");
        }
        sb.append(p.getScaledPhred(q, 0, arm));
      }
      sb.append(StringUtils.LS).append("    ");
    }
    sb.append(StringUtils.LS);
    sb.append(StringUtils.LS);

    if (p.isCG()) {
      sb.append("overlap = ").append(toString(p.overlapDistribution())).append(StringUtils.LS);
      sb.append("gap = ").append(toString(p.gapDistribution())).append(StringUtils.LS);
      sb.append("smallgap = ").append(toString(p.smallGapDistribution())).append(StringUtils.LS);
      sb.append("overlap2 = ").append(toString(p.overlapDistribution2())).append(StringUtils.LS);
    }
    return sb.toString();
  }

  /**
   * Dump string form of genome priors
   * @param p the genome priors
   * @return string of genome priors
   */
  public static String toGenomePriorProperties(GenomePriorParams p) {
    final StringBuilder sb = new StringBuilder();
    sb.append(String.format(Locale.ROOT, "genome_snp_rate_hetero = %.6f" + StringUtils.LS, p.genomeSnpRate(true)));
    sb.append(String.format(Locale.ROOT, "genome_snp_rate_homo = %.6f" + StringUtils.LS, p.genomeSnpRate(false)));
    sb.append(StringUtils.LS);

    sb.append(String.format(Locale.ROOT, "genome_mnp_base_rate_hetero = %.6f" + StringUtils.LS, p.genomeMnpBaseRate(true)));
    sb.append(String.format(Locale.ROOT, "genome_mnp_base_rate_homo = %.6f" + StringUtils.LS, p.genomeMnpBaseRate(false)));
    sb.append("genome_mnp_distribution = ").append(toString(p.genomeMnpDistribution())).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    sb.append(String.format(Locale.ROOT, "genome_indel_event_rate = %.6f" + StringUtils.LS, p.genomeIndelEventRate()));
    sb.append(String.format(Locale.ROOT, "genome_indel_event_fraction = %.6f" + StringUtils.LS, p.genomeIndelEventFraction()));
    sb.append(String.format(Locale.ROOT, "genome_indel_length_decay = %.6f" + StringUtils.LS, p.genomeIndelLengthDecay()));
    sb.append("genome_indel_distribution = ").append(toString(p.genomeIndelDistribution())).append(StringUtils.LS);
    sb.append(StringUtils.LS);

    final HypothesesSnp h = new HypothesesSnp(SimplePossibility.SINGLETON, p, false, -1);
    for (int i = 0; i < h.size(); ++i) {
      final String call = h.code().homozygous(i) ? h.description().name(i) : h.name(i);
      final String call2 = call.toLowerCase(Locale.ROOT).replace(VariantUtils.COLON, '_');
      final double[] dist = p.getPriorDistr(call);
      for (int b = 0; b < 4; ++b) {
        sb.append(String.format(Locale.ROOT, "%s_%s = %.8f" + StringUtils.LS, Character.toLowerCase(DnaUtils.getBase(b + 1)), call2, dist[b]));
      }
      sb.append(StringUtils.LS);
    }
    return sb.toString();
  }

  // See test for details of these numbers
  private static final int[] PHRED_FROM_N = {3/* 2 */, 2, 1, 1, 1, 1, 1, 1/* 9 */};

  /**
   * Given the number of mappings for a read return an appropriate phred score for
   * the probability that it is incorrect.
   * @param n number of times read is mapped.
   * @return phred score.
   */
  public static int phredFromN(final int n) {
    if (n <= 1) {
      throw new IllegalArgumentException("Invalid record count:" + n);
    }
    if (n > 9) {
      return 0;
    }
    return PHRED_FROM_N[n - 2];
  }

  /**
   * Return the previous reference nucleotide as a character
   * @param reference byte array of the reference
   * @param position position to return the previous nucleotide of
   * @return previous reference nucleotide character, or N if out of range.
   */
  public static char getPreviousRefNt(byte[] reference, int position) {
    return position <= 0 || position + 1 > reference.length ? 'N' : DnaUtils.getBase(reference[position - 1]);
  }

  /**
   * Return a normalised array of possibilities
   * @param possibilities the possibilities which are to be normalised
   * @param arith the arithmetic
   * @return the normalised possibilities array
   */
  public static double[] normalisePossibilities(double[] possibilities, PossibilityArithmetic arith) {
    final double[] normalPosses = new double[possibilities.length];
    double sum = arith.zero();
    for (final double poss : possibilities) {
      sum = arith.add(sum, poss);
    }
    for (int i = 0; i < possibilities.length; ++i) {
      final double poss = arith.divide(possibilities[i], sum);
      normalPosses[i] = poss;
    }
    return normalPosses;
  }

  /**
   * Convert an array of probabilities into possibilities
   * @param probabilities the probabilities which are to be converted
   * @param arith the arithmetic
   * @return the normalised possibilities array
   */
  public static double[] prob2Poss(double[] probabilities, PossibilityArithmetic arith) {
    final double[] posses = new double[probabilities.length];
    for (int i = 0; i < probabilities.length; ++i) {
      posses[i] = arith.prob2Poss(probabilities[i]);
    }
    return posses;
  }

  /**
   * Take a string representing a variation and split it into two diploid halves (these
   * will be identical if the call is homozygous). The returned halves also have any "D"s
   * or "i"s removed.
   * @param str a variation which can be heterozygous or homozygous and contain "D" or "i".
   * @return a string array containing the left([0]) and right([1]) predictions (it is always of length 2).
   */
  public static String[] normalizePair(final String str) {
    final String[] res = new String[2];
    if (str.indexOf(VariantUtils.COLON) != -1) {
      //heterozygous
      final String[] split = StringUtils.split(str, VariantUtils.COLON);
      if (split.length != 2) {
        throw new IllegalArgumentException("Invalid variation: " + str);
      }
      res[0] = split[0];
      res[1] = split[1];
    } else {
      //homozygous
      res[0] = str;
      res[1] = str;
    }
    assert res[0] != null && res[1] != null;
    return res;
  }
}
