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
package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.ReferenceBasedFactory;

/**
 * Selects between one factory or another factory depending on the position being populated.
 */
public class SwitchingReferenceBasedBuffer<D> extends ReferenceBasedBuffer<D> {

  private final ReferenceBasedFactory<D> mSecondFactory;
  private final int mSwitchPoint;

  /**
   * @param firstFactory for creating objects in buffer.
   * @param secondFactory for creating objects in buffer.
   * @param switchPoint the first coordinate to get objects made by the second factory.
   * @param template nucleotides.
   * @param start of region being processed on template (0 based)
   */
  public SwitchingReferenceBasedBuffer(ReferenceBasedFactory<D> firstFactory, ReferenceBasedFactory<D> secondFactory, int switchPoint, byte[] template, int start) {
    super(firstFactory, template, start);    //To change body of overridden methods use File | Settings | File Templates.
    mSecondFactory = secondFactory;
    mSwitchPoint = switchPoint;
  }

  @Override
  protected D make(int index) {
    final int nt = mTemplate[index] - 1;
    assert -1 <= nt && nt < 4;
    return index < mSwitchPoint ? mFactory.make(nt) : mSecondFactory.make(nt);
  }
}
