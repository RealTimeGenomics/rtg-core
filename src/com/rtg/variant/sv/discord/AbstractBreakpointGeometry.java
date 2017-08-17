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

import com.rtg.util.StringUtils;
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
  abstract int getXLo(); // AKA x()
  abstract int getXHi(); // AKA z()

  abstract String getYName();
  abstract int getYLo(); // AKA y()
  abstract int getYHi(); // AKA w()

  protected abstract int getRLo(); // AKA r()
  protected abstract int getRHi(); // AKA s()

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

    if (!over(this.getXLo(), this.getXHi(), that.getXLo(), that.getXHi())) {
      return false;
    }

    if (!over(this.getYLo(), this.getYHi(), that.getYLo(), that.getYHi())) {
      return false;
    }

    if (!over(this.getRLo(), this.getRHi(), that.getRLo(), that.getRHi())) {
      return false;
    }

    return true;
  }

  static Interval intersect(int a0, int b0, int a1, int b1) {
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

    final Interval x = intersect(this.getXLo(), this.getXHi(), that.getXLo(), that.getXHi());
    if (x == null) {
      return null;
    }
    final Interval y = intersect(this.getYLo(), this.getYHi(), that.getYLo(), that.getYHi());
    if (y == null) {
      return null;
    }
    final Interval r = intersect(this.getRLo(), this.getRHi(), that.getRLo(), that.getRHi());
    if (r == null) {
      return null;
    }

    //TODO prove that xp, yp, rp cannot be null
    final Interval xp = intersect(x.getA(), x.getB(), x(r.getA(), y.getB()), x(r.getB(), y.getA()));
    final Interval yp = intersect(y.getA(), y.getB(), y(r.getA(), x.getB()), y(r.getB(), x.getA()));
    final Interval rp = intersect(r.getA(), r.getB(), r(x.getA(), y.getA()), r(x.getB(), y.getB()));
    final AbstractBreakpointGeometry res = makeBreakpointGeometry(xp, yp, rp);
    //System.err.println(" res=" + res);
    assert res.integrity();
    return res;
  }

  private AbstractBreakpointGeometry makeBreakpointGeometry(final Interval xp, final Interval yp, final Interval rp) {
    return new BreakpointGeometry(getOrientation(), getXName(), getYName(), xp.getA(), xp.getB(), yp.getA(), yp.getB(), rp.getA(), rp.getB());
  }

  int r(int x, int y) {
    return getOrientation().r(x, y);
  }

  int x(int r, int y) {
    return getOrientation().x(r, y);
  }

  int y(int r, int x) {
    return getOrientation().y(r, x);
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

    final Interval x = union(this.getXLo(), this.getXHi(), that.getXLo(), that.getXHi());
    final Interval y = union(this.getYLo(), this.getYHi(), that.getYLo(), that.getYHi());
    final Interval r = union(this.getRLo(), this.getRHi(), that.getRLo(), that.getRHi());

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
    sb.append(" x=").append(getXLo()).append(",").append(getXHi()).append(":").append(getXName());
    sb.append(" y=").append(getYLo()).append(",").append(getYHi()).append(":").append(getYName());
    sb.append(" r=").append(getRLo()).append(",").append(getRHi());
  }

  private boolean equalsLocal(Object obj) {
    final AbstractBreakpointGeometry that = (AbstractBreakpointGeometry) obj;
    return this.getXName().equals(that.getXName()) && this.getYName().equals(that.getYName())
        && this.getXLo() == that.getXLo() && this.getYLo() == that.getYLo()
        && this.getXHi() == that.getXHi() && this.getYHi() == that.getYHi()
        && this.getRLo() == that.getRLo() && this.getRHi() == that.getRHi();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof AbstractBreakpointGeometry)) {
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
    final int x = bg.getXLo();
    final int y = bg.getYLo();
    final int r = bg.getRLo();
    final int s = bg.getRHi();
    final Orientation or = bg.getOrientation();
    final int rx = or.x(r - or.y(y));
    final int ry = or.y(r - or.x(x));
    final int w = bg.getYHi();
    final int sx = or.x(s - or.y(w));
    final int z = bg.getXHi();
    final int sy = or.y(s - or.x(z));
    sb.append(gnu(max(or.xDir(), x, rx), y));
    sb.append(gnu(x, max(or.yDir(), y, ry)));
    sb.append(gnu(x, w));
    sb.append(gnu(max(or.xDir(), x, sx), w));
    sb.append(gnu(z, max(or.yDir(), y, sy)));
    sb.append(gnu(z, y));
    sb.append(gnu(max(or.xDir(), x, rx), y));
    if (max(or.xDir(), x, rx) != rx) {
      sb.append(StringUtils.LS);
      sb.append(gnu(rx, y));
      sb.append(gnu(max(or.xDir(), x, rx), y));
    }
    if (max(or.yDir(), y, ry) != ry) {
      sb.append(StringUtils.LS);
      sb.append(gnu(x, ry));
      sb.append(gnu(x, max(or.yDir(), y, ry)));
    }
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
    final int xo = getOrientation().xDir();
    final int yo = getOrientation().yDir();
    for (int xi = getXLo(); isInRange(getXLo(), getXHi(), xi); xi += xo) {
      for (int yi = getYLo(); isInRange(getYLo(), getYHi(), yi); yi += yo) {
        count += probe(xi, yi);
      }
    }
    return count;
  }

  int probe(int x, int y) {
    final int r = r(x, y);
    final boolean pr = isInRange(getXLo(), getXHi(), x) && isInRange(getYLo(), getYHi(), y) && isInRange(getRLo(), getRHi(), r);
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
   * Find the earliest position along the specified axis that is covered by this geometry.
   * @param name of the axis along which a position is being found.
   * @return the position (null if neither axis matches the name).
   */
  Integer position(final String name) {
    if (name.equals(getXName())) {
      return getXLo();
    } else if (name.equals(getYName())) {
      return getYLo();
    } else {
      return null;
    }
  }

  @Override
  public boolean integrity() {
    final Orientation or = getOrientation();
    if (or.xDir() == +1) {
      Exam.assertTrue(getXLo() < getXHi());
    } else {
      Exam.assertTrue(getXLo() > getXHi());
    }
    if (or.yDir() == +1) {
      Exam.assertTrue(getYLo() < getYHi());
    } else {
      Exam.assertTrue(getYLo() > getYHi());
    }
    Exam.assertTrue(this.toString(), getRLo() < getRHi());

    boolean error = false;
    final StringBuilder sb = new StringBuilder();
    //constraints on the range of r
    final int rxy = r(getXLo(), getYLo());
    //Assert.assertTrue(getR() >= rxy);

    final int rxw = r(getXLo(), getYHi());
    final int rzy = r(getXHi(), getYLo());
    final int rmax = Math.min(rxw, rzy);
    //Assert.assertTrue(getR() <= rmax);
    if (!(rxy <= getRLo() && getRLo() <= rmax)) {
      sb.append("RLo=").append(getRLo()).append(" not in range ").append(rxy).append(":").append(rmax).append(LS);
      error = true;
    }

    //constraints on the range of s
    final int rzw = r(getXHi(), getYHi());
    final int smin = Math.max(rxw, rzy);
    if (!(smin <= getRHi() && getRHi() <= rzw)) {
      sb.append("RHi=").append(getRHi()).append(" not in range ").append(smin).append(":").append(rzw).append(LS);
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
