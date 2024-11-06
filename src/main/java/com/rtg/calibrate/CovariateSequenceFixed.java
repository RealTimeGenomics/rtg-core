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

package com.rtg.calibrate;

import java.util.List;
import java.util.Locale;

import htsjdk.samtools.SAMRecord;

/**
 * This covariate grows in size as more sequence names are encountered.
 * (Calling parse or value with a new sequence name adds that sequence name
 * into the set automatically).
 *
 */
public final class CovariateSequenceFixed extends CovariateImpl {

  private final List<String> mSequenceNames;

  /**
   * Construct sequence covariate on a known fixed set of sequence names. It is assumed
   * that the supplied list of names is ordered such that <code>sam.getReferenceIndex()</code>
   * will correctly index <code>sequenceNames</code>.
   * @param sequenceNames the sequence names
   */
  public CovariateSequenceFixed(final List<String> sequenceNames) {
    super(CovariateEnum.SEQUENCE.name().toLowerCase(Locale.ROOT), sequenceNames.size());
    mSequenceNames = sequenceNames;
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    return sam.getReferenceIndex();
  }

  @Override
  public int parse(final String name) {
    final int pos = mSequenceNames.indexOf(name);
    if (pos >= 0) {
      return pos;
    }
    throw new UnsupportedOperationException(name + " cf. " + mSequenceNames);
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.SEQUENCE;
  }

  @Override
  public String valueString(final int v) {
    return mSequenceNames.get(v);
  }

}
