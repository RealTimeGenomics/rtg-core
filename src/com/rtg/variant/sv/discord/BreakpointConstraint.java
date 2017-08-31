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
    this(rec, mo, rgs, 0);
  }

  /**
   * @param rec SAM record for one end of a paired read.
   * @param mo orientation expected by the machine technology.
   * @param rgs statistics for gap length etc.
   * @param overlapFraction fraction of read length to assume can overlap the break point
   */
  BreakpointConstraint(SAMRecord rec, MachineOrientation mo, ReadGroupStats rgs, double overlapFraction) {
    //this(makeGeometry(rec, mo, rgs), rgs.gapMean(), rgs.gapStdDev());
    mProxy = makeGeometry(rec, mo, rgs, overlapFraction);
    mMeanR = r(mProxy.getXLo(), mProxy.getYLo()) + rgs.gapMean();
    mStdDeviation = rgs.fragmentStdDev();

    assert getRLo() < mMeanR && mMeanR < getRHi() : this.toString(); // XXX if read overlaps are very common this may not be true, since RLo is clipped
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

  // Get probable mate alignment end (0-based exclusive)
  static int mateEnd(SAMRecord rec) {
    final Integer me = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_MATE_END);
    if (me == null) {
      // Default assume the mate spans as much template as the read itself does (i.e. similar length and alignment characteristics)
      final int length = rec.getAlignmentEnd() - (rec.getAlignmentStart() - 1);
      return (rec.getMateAlignmentStart() - 1) + length;
    }
    return me;
  }


  private static int overlap(SAMRecord rec, int readLength, double overlapFraction) {
    final Integer ii = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES);
    if (ii == null) {
      //default case when attributes not supplied
      return fixedOverlap(readLength, overlapFraction);
    }
    return ii;
  }

  private static int fixedOverlap(int readLength, double overlapFraction) {
    return (int) (readLength * overlapFraction);
  }

  static BreakpointGeometry makeGeometry(SAMRecord rec, MachineOrientation mo, ReadGroupStats rgs, double overlapFraction) {
    final Orientation orientation = Orientation.orientation(rec, mo);
    final String xName = rec.getReferenceName();
    final String yName = rec.getMateReferenceName();
    final int start = rec.getAlignmentStart() - 1; // To 0-based
    final int end = rec.getAlignmentEnd(); // 1-based inclusive == 0-based exclusive
    final int mateStart = rec.getMateAlignmentStart() - 1;
    final int mateEnd = mateEnd(rec);
    final int breakpointOverlap = overlap(rec, end - start, overlapFraction);
    final int mateBreakpointOverlap = fixedOverlap(mateEnd - mateStart, overlapFraction);
    final int x = orientation.xDir() == +1 ?  end - breakpointOverlap : start + breakpointOverlap;
    final int y = orientation.yDir() == +1 ?  mateEnd - mateBreakpointOverlap : mateStart + mateBreakpointOverlap;
    final int min = rgs.gapMin();
    final int max = rgs.gapMax();
    assert min < rgs.gapMean() && rgs.gapMean() < max;
    final int xHi = x + orientation.x(max);
    final int yHi = y + orientation.y(max);
    final int rLo = orientation.r(x, y) + Math.max(min, 0); // rLo can't be left of x,y
    final int rHi = orientation.r(x, y) + max;
    return new BreakpointGeometry(orientation, xName, yName, x, xHi, y, yHi, rLo, rHi);
  }

  /**
   * Check if this is a constraint that fits within the ranges for a normal
   * non-discordant case.
   * @param rgs read group statistics
   * @return true iff not discordant;
   */
  boolean isConcordant(ReadGroupStats rgs) {
    return getXName().equals(getYName())
      && (getOrientation() == Orientation.UD || getOrientation() == Orientation.DU)
      && (getOrientation().r(getXLo(), getYLo()) + rgs.gapMin()) <= 0
      && (getOrientation().r(getXLo(), getYLo()) + rgs.gapMax()) >= 0;
  }

  /**
   * Get an object representing the break point position.
   * @return the position
   */
  BreakpointPosition position() {
    final int rm = (int) MathUtils.round(mMeanR);
    final int r = rm >= getRHi() ? getRHi() - 1 : rm <= getRLo() ? getRLo() : rm;
    //System.err.println("rm=" + rm + " r=" + r + " ox=" + ox + " xy=" + xy + " z=" + z);

    final int ax = min(getOrientation().xDir(), x(r, getYLo()), getXHi());
    final int bx = min(-getOrientation().xDir(), x(r, getYHi()), getXLo());
    final int mid = (ax + bx) / 2;

    final int ay = min(getOrientation().yDir(), y(r, getXLo()), getYHi());
    final int by = min(-getOrientation().yDir(), y(r, getXHi()), getYLo());
    final int midy = (ay + by) / 2;

    //System.err.println("a=" + a + " b=" + b);
    return (ax < bx) ? new BreakpointPosition(ax, mid, bx, midy) : new BreakpointPosition(bx, mid, ax, midy);
  }

  private static int min(int or, int a, int b) {
    assert or == +1 || or == -1;
    if (or == +1) {
      return Math.min(a, b);
    }
    return Math.max(a, b);
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

  String gnuPlotMean() {
    final double r = mMeanR;
    return gnu(getXLo(), y(r, getXLo())) + gnu(x(r, getYLo()), getYLo());
  }

  private double x(double r, int y) {
    final Orientation or = getOrientation();
    return or.x(r - or.y(y));
  }

  private double y(double r, int x) {
    final Orientation or = getOrientation();
    return or.y(r - or.x(x));
  }

  /**
   * @param that constraint being added.
   * @param geometry the purely geometric part which has already been constructed.
   * @return the fully fleshed out <code>BreakpointConstraint</code>.
   */
  private BreakpointConstraint makeDistribution(BreakpointConstraint that, final AbstractBreakpointGeometry geometry) {
    //see discordTheory documentation
    final double w1 = 0.5 / (rStdDev() * rStdDev());
    assert w1 > 0.0 && Double.isFinite(w1);
    final double w2 = 0.5 / (that.rStdDev() * that.rStdDev());
    assert w2 > 0.0 && Double.isFinite(w2);
    final double w = w1 + w2;
    assert w > 0.0 && Double.isFinite(w);
    final double stdDev = Math.sqrt(0.5 / w);
    assert stdDev > 0.0 && Double.isFinite(stdDev);
    final double v1 = rMean() * w1;
    assert Double.isFinite(v1) : v1;
    final double v2 = that.rMean() * w2;
    assert Double.isFinite(v2);
    final double v = v1 + v2;
    assert Double.isFinite(v);
    final double mean = v / w;
    assert Double.isFinite(mean) : mean;
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

  @Override
  BreakpointConstraint flip() {
    return new BreakpointConstraint(mProxy.flip(), mMeanR, mStdDeviation);
  }

  @Override
  protected Orientation getOrientation() {
    return mProxy.getOrientation();
  }

  @Override
  protected int getRLo() {
    return mProxy.getRLo();
  }

  @Override
  protected int getRHi() {
    return mProxy.getRHi();
  }

  @Override
  protected int getYHi() {
    return mProxy.getYHi();
  }

  @Override
  protected int getXLo() {
    return mProxy.getXLo();
  }

  @Override
  protected String getXName() {
    return mProxy.getXName();
  }

  @Override
  protected int getYLo() {
    return mProxy.getYLo();
  }

  @Override
  protected String getYName() {
    return mProxy.getYName();
  }

  @Override
  protected int getXHi() {
    return mProxy.getXHi();
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
    Exam.assertTrue(Double.isFinite(mMeanR));
    Exam.assertTrue(0.0 < mStdDeviation && Double.isFinite(mStdDeviation));
    return true;
  }
}
