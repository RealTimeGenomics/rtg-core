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

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Describes a six-sided region in two dimensional space.
 * The axes are positions along two (possibly different) sequences being the
 * putative positions of a breakpoint. The six sides are two parallel to the x and y axes and
 * two parallel to a diagonal (which diagonal depends on the orientation).
 * <br>
 * This class is used to represent individual break-point constraints determined by
 * paired-end reads and unions and intersections of these as well as geometrical transformations
 * of these.
 * <br>
 * An illustration of the conventions used in representing these shapes is given below:
 * <img src="doc-files/BreakpointDiagrams/Slide6.jpg" alt="image">
 */
public abstract class AbstractBreakpointGeometry extends IntegralAbstract {

  /**
   * Check if a one-dimensional interval overlaps. a0, b0 and a1, b1 form
   * pairs that can be oriented positively (a &lt; b) or negatively (a &gt; b).
   * a = b implies an empty interval (because the b's are exclusive).
   * @param a0 start co-ordinate for first interval (inclusive).
   * @param b0 end co-ordinate for first interval(exclusive).
   * @param a1 start co-ordinate for second interval (inclusive).
   * @param b1 length of second interval (0 based, exclusive).
   * @return true iff the the two intervals overlap.
   */
  protected static boolean over(int a0, int b0, int a1, int b1) {
    if (b0 > a0) {
      return b1 > a0 && b0 > a1 && a1 < b1;
    } else if (a0 > b0) {
      return b1 < a0 && b0 < a1 && a1 > b1;
    } else {
      return false;
    }
  }

  abstract AbstractBreakpointGeometry flip();

  abstract Orientation getOrientation();

  abstract String getXName();
  abstract int getX();
  abstract int getZ();

  abstract String getYName();
  abstract int getY();
  abstract int getW();

  protected abstract int getR();
  protected abstract int getS();

  protected boolean overlap(AbstractBreakpointGeometry that) {
    if (this.getOrientation() != that.getOrientation()) {
      return false;
    }
    if (!this.getXName().equals(that.getXName())) {
      return false;
    }
    if (!this.getYName().equals(that.getYName())) {
      return false;
    }

    if (!over(this.getX(), this.getZ(), that.getX(), that.getZ())) {
      return false;
    }

    if (!over(this.getY(), this.getW(), that.getY(), that.getW())) {
      return false;
    }

    if (!over(this.getR(), this.getS(), that.getR(), that.getS())) {
      return false;
    }

    return true;
  }

  static Interval inter(int a0, int b0, int a1, int b1) {
    //System.err.println("inter a0=" + a0 + " b0=" + b0 + " a1=" + a1 + " b1=" + b1);
    final int a2;
    final int b2;
    if (a0 < b0) {
      if (a1 >= b1) {
        return null;
      }
      a2 = Math.max(a0, a1);
      b2 = Math.min(b0, b1);
      if (a2 >= b2) {
        return null;
      }
    } else if (a0 > b0) {
      if (a1 <= b1) {
        return null;
      }
      a2 = Math.min(a0, a1);
      b2 = Math.max(b0, b1);
      if (a2 <= b2) {
        return null;
      }
    } else {
      return null;
    }
    assert a2 != b2;
    //System.err.println("inter res=" + res);
    return new Interval(a2, b2);
  }

  AbstractBreakpointGeometry intersect(AbstractBreakpointGeometry that) {
    //System.err.println("this=" + this);
    //System.err.println("that=" + that);
    if (this.getOrientation() != that.getOrientation()) {
      return null;
    }
    if (!this.getXName().equals(that.getXName())) {
      return null;
    }
    if (!this.getYName().equals(that.getYName())) {
      return null;
    }

    final Interval x = inter(this.getX(), this.getZ(), that.getX(), that.getZ());
    if (x == null) {
      return null;
    }
    final Interval y = inter(this.getY(), this.getW(), that.getY(), that.getW());
    if (y == null) {
      return null;
    }
    final Interval r = inter(this.getR(), this.getS(), that.getR(), that.getS());
    if (r == null) {
      return null;
    }

    //TODO prove that xp, yp, rp cannot be null
    final Interval xp = inter(x.getA(), x.getB(), x(r.getA(), y.getB()), x(r.getB(), y.getA()));
    final Interval yp = inter(y.getA(), y.getB(), y(r.getA(), x.getB()), y(r.getB(), x.getA()));
    final Interval rp = inter(r.getA(), r.getB(), r(x.getA(), y.getA()), r(x.getB(), y.getB()));

    final AbstractBreakpointGeometry res = makeBreakpointGeometry(xp, yp, rp);
    //System.err.println(" res=" + res);
    assert res.integrity();
    return res;
  }

  private AbstractBreakpointGeometry makeBreakpointGeometry(final Interval xp, final Interval yp, final Interval rp) {
    return new BreakpointGeometry(getOrientation(), getXName(), getYName(), xp.getA(), xp.getB(), yp.getA(), yp.getB(), rp.getA(), rp.getB());
  }

  int r(int x, int y) {
    final Orientation or = getOrientation();
    return or.x(x) + or.y(y);
  }

  int x(int r, int y) {
    final Orientation or = getOrientation();
    return or.x(r - or.y(y));
  }

  int y(int r, int x) {
    final Orientation or = getOrientation();
    return or.y(r - or.x(x));
  }

  AbstractBreakpointGeometry union(AbstractBreakpointGeometry that) {
    if (this.getOrientation() != that.getOrientation()) {
      return null;
    }
    if (!this.getXName().equals(that.getXName())) {
      return null;
    }
    if (!this.getYName().equals(that.getYName())) {
      return null;
    }

    final Interval x = union(this.getX(), this.getZ(), that.getX(), that.getZ());
    final Interval y = union(this.getY(), this.getW(), that.getY(), that.getW());
    final Interval r = union(this.getR(), this.getS(), that.getR(), that.getS());

    final AbstractBreakpointGeometry res = makeBreakpointGeometry(x, y, r);
    //System.err.println(this.gnuPlot());
    //System.err.println(that.gnuPlot());
    //System.err.println(res.gnuPlot());

    //System.err.println("res=" + res);
    assert res.integrity();
    return res;
  }

  static Interval union(int a0, int b0, int a1, int b1) {
    final int a2;
    final int b2;
    if (a0 < b0) {
      assert a1 < b1;
      a2 = Math.min(a0, a1);
      b2 = Math.max(b0, b1);
      assert a2 < b2;
    } else if (a0 > b0) {
      assert a1 > b1;
      a2 = Math.max(a0, a1);
      b2 = Math.min(b0, b1);
      assert a2 > b2;
    } else {
      throw new RuntimeException();
    }
    return new Interval(a2, b2);
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Break-point constraint:").append(getOrientation());
    sb.append(" x=").append(getX()).append(",").append(getZ()).append(":").append(getXName());
    sb.append(" y=").append(getY()).append(",").append(getW()).append(":").append(getYName());
    sb.append(" r=").append(getR()).append(",").append(getS());
  }

  private boolean equalsLocal(Object obj) {
    final AbstractBreakpointGeometry that = (AbstractBreakpointGeometry) obj;
    return this.getXName().equals(that.getXName()) && this.getYName().equals(that.getYName())
        && this.getX() == that.getX() && this.getY() == that.getY()
        && this.getZ() == that.getZ() && this.getW() == that.getW()
        && this.getR() == that.getR() && this.getS() == that.getS();
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    return this.equalsLocal(obj) || flip().equalsLocal(obj);
  }

  @Override
  public int hashCode() {
    // pacify findbugs
    return super.hashCode();
  }

  /**
   * Generate a graph of the periphery of the geometry suitable for input
   * to gnuplot (that is x,y pairs of the corner points).
   * @return a string for input to gnuplot.
   */
  public String gnuPlot() {
    return gnuPlot(this);
  }

  static int max(final int or, final int a, final int b) {
    if (or == +1) {
      return Math.max(a, b);
    }
    assert or == -1;
    return Math.min(a, b);
  }

  /**
   * Generate a graph of the periphery of the geometry suitable for input
   * to gnuplot (that is x,y pairs of the corner points).
   * @param bg the geometry to plot
   * @return a string for input to gnuplot.
   */
  static String gnuPlot(final AbstractBreakpointGeometry bg) {
    final StringBuilder sb = new StringBuilder();
    final int x = bg.getX();
    final int y = bg.getY();
    final int r = bg.getR();
    final int s = bg.getS();
    final Orientation or = bg.getOrientation();
    final int rx = or.x(r - or.y(y));
    final int ry = or.y(r - or.x(x));
    final int w = bg.getW();
    final int sx = or.x(s - or.y(w));
    final int z = bg.getZ();
    final int sy = or.y(s - or.x(z));
    sb.append(gnu(max(or.getX(), x, rx), y));
    sb.append(gnu(x, max(or.getY(), y, ry)));
    sb.append(gnu(x, w));
    sb.append(gnu(max(or.getX(), x, sx), w));
    sb.append(gnu(z, max(or.getY(), y, sy)));
    sb.append(gnu(z, y));
    sb.append(gnu(max(or.getX(), x, rx), y));
    return sb.toString();
  }

  /**
   * One line of gnuplot output.
   * @param x co-ordinate of point.
   * @param y co-ordinate of point.
   * @return one line of gnuplot output.
   */
  static String gnu(double x, double y) {
    return Utils.realFormat(x, 1) + "\t " + Utils.realFormat(y, 1) + LS;
  }

  @Override
  public boolean globalIntegrity() {
    //System.err.println(this);
    final int count = count();
    //System.err.println("count=" + count);
    integrity();
    Exam.assertTrue(count > 0);
    return true;
  }

  /**
   * Count the number of points in the region.
   * @return the count.
   */
  int count() {
    int count = 0;
    final int xo = getOrientation().getX();
    final int yo = getOrientation().getY();
    for (int xi = getX(); isInRange(getX(), getZ(), xi); xi += xo) {
      for (int yi = getY(); isInRange(getY(), getW(), yi); yi += yo) {
        count += probe(xi, yi);
      }
    }
    return count;
  }

  int probe(int x, int y) {
    final int r = r(x, y);
    final boolean pr = isInRange(getX(), getZ(), x) && isInRange(getY(), getW(), y) && isInRange(getR(), getS(), r);
    return pr ? 1 : 0;
  }

  /**
   * Check that x is in the range from a (inclusive) to b (exclusive).
   * Works for either a &lt; b or a &gt; b.
   * @param a start of range (inclusive).
   * @param b end of range (exclusive).
   * @param x point to be checked
   * @return true iff x is range from a to b.
   */
  static boolean isInRange(final int a, final int b, final int x) {
    //System.err.println("probe a=" + a + " b=" + b + " x=" + x);
    if (a < b) {
      return a <= x && x < b;
    } else {
      assert b < a;
      return a >= x && x > b;
    }
  }

  /**
   * Check that x which is assumed to be an exclusive bound is in the range from a (inclusive) to b (exclusive bound but inclusive because checking another bound).
   * Works for either a &lt; b or a &gt; b.
   * @param a start of range (inclusive).
   * @param b end of range (exclusive).
   * @param x exclusive bound to be checked
   * @return true iff x is range from a to b.
   */
  static boolean isInRangeLimit(final int a, final int b, final int x) {
    //System.err.println("probe a=" + a + " b=" + b + " x=" + x);
    if (a < b) {
      return a <= x && x <= b;
    } else {
      assert b < a;
      return a >= x && x >= b;
    }
  }

  /**
   * Find the earliest position along the specified axis that is covered by this geometry.
   * @param name of the axis along which a position is being found.
   * @return the position (null if neither axis matches the name).
   */
  Integer position(final String name) {
    if (name.equals(getXName())) {
      return getX();
    } else if (name.equals(getYName())) {
      return getY();
    } else {
      return null;
    }
  }

  @Override
  public boolean integrity() {
    final Orientation or = getOrientation();
    if (or.getX() == +1) {
      Exam.assertTrue(getX() < getZ());
    } else {
      Exam.assertTrue(getX() > getZ());
    }
    if (or.getY() == +1) {
      Exam.assertTrue(getY() < getW());
    } else {
      Exam.assertTrue(getY() > getW());
    }
    Exam.assertTrue(this.toString(), getR() < getS());

    boolean error = false;
    final StringBuilder sb = new StringBuilder();
    //constraints on the range of r
    final int rxy = r(getX(), getY());
    //Assert.assertTrue(getR() >= rxy);

    final int rxw = r(getX(), getW());
    final int rzy = r(getZ(), getY());
    final int rmax = Math.min(rxw, rzy);
    //Assert.assertTrue(getR() <= rmax);
    if (!(rxy <= getR() && getR() <= rmax)) {
      sb.append("R=").append(getR()).append(" not in range ").append(rxy).append(":").append(rmax).append(LS);
      error = true;
    }

    //constraints on the range of s
    final int rzw = r(getZ(), getW());
    final int smin = Math.max(rxw, rzy);
    if (!(smin <= getS() && getS() <= rzw)) {
      sb.append("S=").append(getS()).append(" not in range ").append(smin).append(":").append(rzw).append(LS);
      error = true;
    }
    //Assert.assertTrue(smin <= getS());
    //Assert.assertTrue(getS() <= rzw);

    //Assert.assertTrue("" + rxw + ":" + r(getZ(), getW()) + ":" + getS(), isInRangeLimit(rxw, r(getZ(), getW()), getS()));
    //Assert.assertTrue(isInRangeLimit(rzy, r(getZ(), getW()), getS()));
    if (error) {
      System.err.println(sb.toString());
      Exam.assertTrue(false);
    }
    return true;
  }



}
