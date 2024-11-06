/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
   * @param initialCapacity starting buffer length
   * @param firstFactory for creating objects in buffer.
   * @param secondFactory for creating objects in buffer.
   * @param switchPoint the first coordinate to get objects made by the second factory.
   * @param template nucleotides.
   * @param start of region being processed on template (0 based)
   */
  public SwitchingReferenceBasedBuffer(int initialCapacity, ReferenceBasedFactory<D> firstFactory, ReferenceBasedFactory<D> secondFactory, int switchPoint, byte[] template, int start) {
    super(initialCapacity, firstFactory, template, start);    //To change body of overridden methods use File | Settings | File Templates.
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
