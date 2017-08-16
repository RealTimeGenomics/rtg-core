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
public final class FlippedProxyBreakpointGeometry extends AbstractBreakpointGeometry {

  private final AbstractBreakpointGeometry mProxy;

  /**
   * @param br the proxy.
   */
  protected FlippedProxyBreakpointGeometry(AbstractBreakpointGeometry br) {
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
  protected int getXLo() {
    return mProxy.getYLo();
  }

  @Override
  protected String getXName() {
    return mProxy.getYName();
  }

  @Override
  protected int getYLo() {
    return mProxy.getXLo();
  }

  @Override
  protected String getYName() {
    return mProxy.getXName();
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
    return mProxy.getXHi();
  }

  @Override
  protected int getXHi() {
    return mProxy.getYHi();
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
