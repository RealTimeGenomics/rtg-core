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

import com.rtg.util.integrity.Exam;

/**
 * Describes a six-sided region in two dimensional space with six explicit values.
 *
 * This class is used to represent individual break-point constraints determined by
 * paired-end reads and unions and intersections of these.
 */
public final class BreakpointGeometry extends AbstractBreakpointGeometry {

  /**
   * Create a break point given the real co-ordinates.
   * Useful for re-creating data from text files etc.
   * @param xName name for x co-ordinate.
   * @param x co-ordinate.
   * @param z co-ordinate.
   * @param yName name for y co-ordinate.
   * @param y co-ordinate.
   * @param w co-ordinate.
   * @param r co-ordinate.
   * @param s co-ordinate.
   * @return the break point.
   */
  public static BreakpointGeometry makeGeometry(final String xName, final int x, final int z, final String yName, final int y, final int w, final int r, final int s) {
    final Orientation orientation = Orientation.orientation(z - x, w - y);
    final BreakpointGeometry bg = new BreakpointGeometry(orientation, xName, yName, x, z, y, w, r, s);
    bg.globalIntegrity();
    return bg;
  }

  private final Orientation mOrientation;
  private final String mXName;
  private final String mYName;
  private final int mX;
  private final int mZ;
  private final int mY;
  private final int mW;
  private final int mR;
  private final int mS;

  private final AbstractBreakpointGeometry mFlip;

  /**
   * @param orientation along the two axes.
   * @param xName name of x sequence.
   * @param yName name of y sequence.
   * @param x first x co-ordinate value
   * @param z second x co-ordinate value
   * @param y first y co-ordinate value
   * @param w second y co-ordinate value
   * @param r first diagonal value
   * @param s second diagonal value
   */
  protected BreakpointGeometry(Orientation orientation, String xName, String yName, int x, int z, int y, int w, int r, int s) {
    super();
    mOrientation = orientation;
    mXName = xName;
    mYName = yName;
    mX = x;
    mZ = z;
    mY = y;
    mW = w;
    mR = r;
    mS = s;
    mFlip = new FlippedProxyBreakpointConstraint(this);
    //assert globalIntegrity();
  }

  /**
   * Get orientation.
   * @return Returns the orientation.
   */
  @Override
  protected Orientation getOrientation() {
    return mOrientation;
  }

  /**
   * Get xName.
   * @return Returns the xName.
   */
  @Override
  protected String getXName() {
    return mXName;
  }

  /**
   * Get yName.
   * @return Returns the yName.
   */
  @Override
  protected String getYName() {
    return mYName;
  }

  /**
   * Get x.
   * @return Returns the x.
   */
  @Override
  protected int getX() {
    return mX;
  }

  /**
   * Get z.
   * @return Returns the z.
   */
  @Override
  protected int getZ() {
    return mZ;
  }

  /**
   * Get y.
   * @return Returns the y.
   */
  @Override
  protected int getY() {
    return mY;
  }

  /**
   * Get w.
   * @return Returns the w.
   */
  @Override
  protected int getW() {
    return mW;
  }

  /**
   * Get r.
   * @return Returns the r.
   */
  @Override
  protected int getR() {
    return mR;
  }

  /**
   * Get s.
   * @return Returns the s.
   */
  @Override
  protected int getS() {
    return mS;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertNotNull(mOrientation);
    Exam.assertNotNull(mXName);
    Exam.assertTrue(0 <= mX);
    Exam.assertNotNull(mYName);
    Exam.assertTrue(0 <= mY);
    return true;
  }

  @Override
  AbstractBreakpointGeometry flip() {
    return mFlip;
  }
}
