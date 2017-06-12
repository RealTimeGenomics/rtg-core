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
package com.rtg.variant.sv.discord;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.sam.SamUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.variant.sv.ReadGroupStats;

import htsjdk.samtools.SAMRecord;

/**
 * Specifies the constraints on the positions of a breakpoint given a single paired end read
 * with both ends successfully mapped.
 * <br><br>
 * The diagrams below illustrate the conventions used in representing the possible range of values
 * (via <code>gapMean()</code> and <code>gapStdDEv()</code>) and the various reads alignments that can give rise to the different cases.
 *
 * <img src="doc-files/BreakpointDiagrams/Slide5.jpg" alt="image">
 * <img src="doc-files/BreakpointDiagrams/Slide1.jpg" alt="image">
 * <img src="doc-files/BreakpointDiagrams/Slide2.jpg" alt="image">
 * <img src="doc-files/BreakpointDiagrams/Slide3.jpg" alt="image">
 * <img src="doc-files/BreakpointDiagrams/Slide4.jpg" alt="image">
 * <img src="doc-files/BreakpointDiagrams/Slide7.jpg" alt="image">
 *
 */
@TestClass({"com.rtg.variant.sv.discord.BreakpointConstraintTest", "com.rtg.variant.sv.discord.FlipTest"})
public final class BreakpointConstraint extends AbstractBreakpointGeometry {

  /** The number of standard deviations delineating concordant vs discordant fragment lengths */
  private static final double DISCORDANT_NUM_STD_DEV = GlobalFlags.getDoubleValue(CoreGlobalFlags.SV_DISCORDANT_STD_DEV);

  /**
   * Get the amount of deviation from the mean considered concordant
   * @param rgs read group statistics.
   * @return the concordant deviation.
   */
  static double concordantDeviation(ReadGroupStats rgs) {
    return DISCORDANT_NUM_STD_DEV * rgs.fragmentStdDev();
  }

  /**
   * Get the minimum gap between first and second read.
   * @param rgs read group statistics.
   * @return the minimum gap.
   */
  static int gapMin(ReadGroupStats rgs) {
    return Math.max(0, (int) (rgs.gapMean() - concordantDeviation(rgs) + .5));
  }

  /**
   * Get the maximum gap between first and second read.
   * @param rgs read group statistics.
   * @return the maximum gap.
   */
  static int gapMax(ReadGroupStats rgs) {
    return (int) (rgs.gapMean() + concordantDeviation(rgs) + .5);
  }


  private final AbstractBreakpointGeometry mProxy;

  /** Mean position of breakpoint in R co-ordinate. */
  private final double mMeanR;
  private final double mStdDeviation;

  /**
   * @param rec SAM record for one end of a paired read.
   * @param mo orientation expected by the machine technology.
   * @param rgs statistics for gap length etc.
   */
  BreakpointConstraint(SAMRecord rec, MachineOrientation mo, ReadGroupStats rgs) {
    //this(makeGeometry(rec, mo, rgs), rgs.gapMean(), rgs.gapStdDev());
    mProxy = makeGeometry(rec, mo, rgs);
    mMeanR = r(mProxy.getX(), mProxy.getY()) + rgs.gapMean();
    mStdDeviation = rgs.fragmentStdDev();
    assert getR() < mMeanR && mMeanR < getS() : this.toString();
    assert globalIntegrity();
  }

  //Used internally and in testing
  BreakpointConstraint(AbstractBreakpointGeometry proxy, double rMean, double rStdDev) {
    mProxy = proxy;
    mMeanR = rMean;
    mStdDeviation = rStdDev;
    //System.err.println(this);
    assert globalIntegrity();
  }

  double x(double r, int y) {
    final Orientation or = getOrientation();
    return or.x(r - or.y(y));
  }

  double y(double r, int x) {
    final Orientation or = getOrientation();
    return or.y(r - or.x(x));
  }

  String gnuPlotMean() {
    final double r = mMeanR;
    return gnu(getX(), y(r, getX())) + gnu(x(r, getY()), getY());
  }

  static int mateEnd(SAMRecord rec) {
    final Integer me = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_MATE_END);
    if (me == null) {
      // Default assume the mate spans as much template as the read itself does (i.e. similar length and alignment characteristics)
      final int length = rec.getAlignmentEnd() - rec.getAlignmentStart();
      return rec.getMateAlignmentStart() + length;
    }
    return me;
  }


  private static int overlap(SAMRecord rec, int readLength) {
    final Integer ii = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES);
    if (ii == null) {
      //default case when attributes not supplied
      return fixedOverlap(readLength);
    }
    return ii;
  }

  private static int fixedOverlap(int readLength) {
    return readLength * ReadGroupStats.ALIGNMENT_IGNORED_FRACTION / 100;
  }

  static BreakpointGeometry makeGeometry(SAMRecord rec, MachineOrientation mo, ReadGroupStats rgs) {
    final Orientation orientation = Orientation.orientation(rec, mo);
    final int start = rec.getAlignmentStart();
    final int end = rec.getAlignmentEnd();
    final int mateStart = rec.getMateAlignmentStart();
    final int mateEnd = mateEnd(rec);
    final int breakpointOverlap = overlap(rec, end - start);
    final int mateBreakpointOverlap = fixedOverlap(mateEnd - mateStart);
    final int x = orientation.getX() == +1 ?  end + 1 - breakpointOverlap : start - 1 + breakpointOverlap;
    final String xName = rec.getReferenceName();
    final String yName = rec.getMateReferenceName();
    final int y = orientation.getY() == +1 ?  mateEnd + 1 - mateBreakpointOverlap : mateStart - 1 + mateBreakpointOverlap;
    final int max = gapMax(rgs);
    final int min = gapMin(rgs);
    assert min < rgs.gapMean() && rgs.gapMean() < max;
    final int z = x + orientation.getX() * max;
    final int w = y + orientation.getY() * max;
    final int r = orientation.x(x) + orientation.y(y) + min;
    final int s = orientation.x(x) + orientation.y(y) + max;
    return new BreakpointGeometry(orientation, xName, yName, x, z, y, w, r, s);
  }

  static int min(int or, int a, int b) {
    assert or == +1 || or == -1;
    if (or == +1) {
      return Math.min(a, b);
    }
    return Math.max(a, b);
  }

  /**
   * Get a VCF record for the specified sequence.
   * @return <code>vcfRecord</code> for displaying the selected breakpoint (null if sequence name doesn't agree).
   */
  BreakpointPosition position() {
    final int rm = (int) MathUtils.round(mMeanR);
    final int r = rm >= getS() ? getS() - 1 : rm <= getR() ? getR() : rm;
    final int xy = x(r, getY());
    final int z = getZ();
    final int ox = getOrientation().getX();
    //System.err.println("rm=" + rm + " r=" + r + " ox=" + ox + " xy=" + xy + " z=" + z);

    final int a = min(ox, xy, z);
    final int b = min(-ox, x(r, getW()), getX());
    final int mid = (a + b) / 2;

    final int ay = min(getOrientation().getY(), y(r, getX()), getW());
    final int by = min(-getOrientation().getY(), y(r, z), getY());
    final int midy = (ay + by) / 2;

    //System.err.println("a=" + a + " b=" + b);
    return (a < b) ? new BreakpointPosition(a, mid, b, midy) : new BreakpointPosition(b, mid, a, midy);
  }

  BreakpointConstraint outputGeometry(final String axis) {
    if (axis.equals(getXName())) {
      return this;
    } else if (axis.equals(getYName())) {
      return flip();
    } else {
      return null;
    }
  }

  /**
   * @param that constraint being added.
   * @param geometry the purely geometric part which has already been constructed.
   * @return the fully fleshed out <code>BreakpointConstraint</code>.
   */
  private BreakpointConstraint makeDistribution(BreakpointConstraint that, final AbstractBreakpointGeometry geometry) {
    //see discordTheory documentation
    final double w1 = 0.5 / (rStdDev() * rStdDev());
    assert w1 > 0.0 && !Double.isInfinite(w1) && !Double.isNaN(w1);
    final double w2 = 0.5 / (that.rStdDev() * that.rStdDev());
    assert w2 > 0.0 && !Double.isInfinite(w2) && !Double.isNaN(w2);
    final double w = w1 + w2;
    assert w > 0.0 && !Double.isInfinite(w) && !Double.isNaN(w);
    final double stdDev = Math.sqrt(0.5 / w);
    assert stdDev > 0.0 && !Double.isInfinite(stdDev) && !Double.isNaN(stdDev);
    final double v1 = rMean() * w1;
    assert !Double.isInfinite(v1) && !Double.isNaN(v1) : v1;
    final double v2 = that.rMean() * w2;
    assert !Double.isInfinite(v2) && !Double.isNaN(v2);
    final double v = v1 + v2;
    assert !Double.isInfinite(v) && !Double.isNaN(v);
    final double mean = v / w;
    assert !Double.isInfinite(mean) && !Double.isNaN(mean) : mean;
    //System.err.println("sigma1=" + gapStdDev() + " sigma2=" + bc.gapStdDev());
    //System.err.println("mu1=" + gapMean() + " mu2=" + bc.gapMean());
    //System.err.println("v1=" + v1 + " v2=" + v2);
    //System.err.println("w1=" + w1 + " w2=" + w2);
    //System.err.println("mean=" + mean + " stdDev=" + stdDev);
    //System.err.println();
    return new BreakpointConstraint(geometry, mean, stdDev);
  }

  @Override
  BreakpointConstraint intersect(AbstractBreakpointGeometry that) {
    final AbstractBreakpointGeometry intersect = mProxy.intersect(that);
    if (intersect == null) {
      return null;
    }
    return makeDistribution((BreakpointConstraint) that, intersect);
  }

  @Override
  BreakpointConstraint union(AbstractBreakpointGeometry that) {
    //System.err.println("union");
    //System.err.println("this=" + this);
    //System.err.println("that=" + that);
    final AbstractBreakpointGeometry union = mProxy.union(that);
    if (union == null) {
      return null;
    }
    return makeDistribution((BreakpointConstraint) that, union);
  }

  /**
   * Check if this is a constraint that fits within the ranges for a normal
   * non-discordant case.
   * @return true iff not discordant;
   */
  boolean isConcordant() {
    if (!getXName().equals(getYName())) {
      return false;
    }
    if (getOrientation() == Orientation.UU || getOrientation() == Orientation.DD) {
      return false;
    }
    return getR() <= 1 && getS() >= 1;
  }

  @Override
  BreakpointConstraint flip() {
    return new BreakpointConstraint(mProxy.flip(), mMeanR, mStdDeviation);
  }

  @Override
  protected Orientation getOrientation() {
    return mProxy.getOrientation();
  }

  @Override
  protected int getR() {
    return mProxy.getR();
  }

  @Override
  protected int getS() {
    return mProxy.getS();
  }

  @Override
  protected int getW() {
    return mProxy.getW();
  }

  @Override
  protected int getX() {
    return mProxy.getX();
  }

  @Override
  protected String getXName() {
    return mProxy.getXName();
  }

  @Override
  protected int getY() {
    return mProxy.getY();
  }

  @Override
  protected String getYName() {
    return mProxy.getYName();
  }

  @Override
  protected int getZ() {
    return mProxy.getZ();
  }

  /**
   * Get the mean gap in r coordinates.
   * This can be used for the length of an indel when
   * <code>isConcordant()</code> is true.
   * @return the mean gap.
   */
  public double rMean() {
    return mMeanR;
  }

  /**
   * The estimated standard deviation of the <code>gapMean()</code>.
   * @return the standard deviation.
   */
  public double rStdDev() {
    return mStdDeviation;
  }

  @Override
  public void toString(StringBuilder sb) {
    super.toString(sb);
    sb.append(" Gap:");
    sb.append(" mean=").append(Utils.realFormat(rMean(), 1));
    sb.append(" std.dev.=").append(Utils.realFormat(rStdDev(), 1));
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(!Double.isInfinite(mMeanR) && !Double.isNaN(mMeanR));
    Exam.assertTrue(0.0 < mStdDeviation && !Double.isInfinite(mStdDeviation) && !Double.isNaN(mStdDeviation));
    return true;
  }
}
