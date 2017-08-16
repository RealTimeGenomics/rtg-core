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
  private final int mXLo;
  private final int mXHi;
  private final int mYLo;
  private final int mYHi;
  private final int mRLo;
  private final int mRHi;

  private final AbstractBreakpointGeometry mFlip;

  /**
   * @param orientation along the two axes.
   * @param xName name of x sequence.
   * @param yName name of y sequence.
   * @param xLo first x co-ordinate value
   * @param xHi second x co-ordinate value
   * @param yLo first y co-ordinate value
   * @param yHi second y co-ordinate value
   * @param rLo first diagonal value
   * @param rHi second diagonal value
   */
  protected BreakpointGeometry(Orientation orientation, String xName, String yName, int xLo, int xHi, int yLo, int yHi, int rLo, int rHi) {
    super();
    mOrientation = orientation;
    mXName = xName;
    mYName = yName;
    mXLo = xLo;
    mXHi = xHi;
    mYLo = yLo;
    mYHi = yHi;
    mRLo = rLo;
    mRHi = rHi;
    mFlip = new FlippedProxyBreakpointGeometry(this);
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
  protected int getXLo() {
    return mXLo;
  }

  /**
   * Get z.
   * @return Returns the z.
   */
  @Override
  protected int getXHi() {
    return mXHi;
  }

  /**
   * Get y.
   * @return Returns the y.
   */
  @Override
  protected int getYLo() {
    return mYLo;
  }

  /**
   * Get w.
   * @return Returns the w.
   */
  @Override
  protected int getYHi() {
    return mYHi;
  }

  /**
   * Get r.
   * @return Returns the r.
   */
  @Override
  protected int getRLo() {
    return mRLo;
  }

  /**
   * Get s.
   * @return Returns the s.
   */
  @Override
  protected int getRHi() {
    return mRHi;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertNotNull(mOrientation);
    Exam.assertNotNull(mXName);
    Exam.assertTrue(0 <= mXLo);
    Exam.assertNotNull(mYName);
    Exam.assertTrue(0 <= mYLo);
    return true;
  }

  @Override
  AbstractBreakpointGeometry flip() {
    return mFlip;
  }
}
