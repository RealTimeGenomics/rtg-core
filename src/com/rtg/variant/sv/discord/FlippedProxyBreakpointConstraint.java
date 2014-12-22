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
import com.rtg.util.integrity.Exam;



/**
 * Flip x and y co-ordinates.
 */
@TestClass(value = {"com.rtg.variant.sv.discord.FlipTest"})
public final class FlippedProxyBreakpointConstraint extends AbstractBreakpointGeometry {

  private final AbstractBreakpointGeometry mProxy;

  /**
   * @param br the proxy.
   */
  protected FlippedProxyBreakpointConstraint(AbstractBreakpointGeometry br) {
    mProxy = br;
  }

  @Override
  AbstractBreakpointGeometry flip() {
    return mProxy;
  }

  @Override
  protected Orientation getOrientation() {
    return mProxy.getOrientation().flip();
  }

  @Override
  protected int getX() {
    return mProxy.getY();
  }

  @Override
  protected String getXName() {
    return mProxy.getYName();
  }

  @Override
  protected int getY() {
    return mProxy.getX();
  }

  @Override
  protected String getYName() {
    return mProxy.getXName();
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
    return mProxy.getZ();
  }

  @Override
  protected int getZ() {
    return mProxy.getW();
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertNotNull(mProxy);
    mProxy.integrity();
    Exam.assertEquals(this.count(), mProxy.count());
    return true;
  }
}
