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

package com.rtg.variant;

import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.ReaderRecord;
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.VariantUtils;

/**
 */
public class SamToMatchCigar extends IntegralAbstract implements SamToMatch {

  private final VariantParams mParams;

  private final ReadParserInterface mParser;

  private final int mQdefault;

  private final MachineErrorChooserInterface mChooser;

  /**
   * @param params command line parameters.
   * @param parser parse and expand cigars.
   * @param chooser gets machine errors for alignment records.
   */
  public SamToMatchCigar(final VariantParams params, final ReadParserInterface parser, final MachineErrorChooserInterface chooser) {
    super();
    mParams = params;
    mParser = parser;
    mQdefault = mParams.qDefault();
    mChooser = chooser;
  }

  @Override
  public boolean process(final byte[] templateBytes, final VariantAlignmentRecord rec) {
    final MachineType machineType = mChooser.machineType(rec.getReadGroup(), rec.isReadPaired());
    final int readScore = VariantUtils.readScoreFromAlignmentRecord(rec, mParams);
    if (readScore < 0) {
      // Invalid read score
      return false;
    }
    try {
      // System.err.println("read=" + read + " cigar=" + cigarString);
      mParser.toMatcher(rec, machineType, mQdefault, templateBytes);
      return true;
    } catch (final BadSuperCigarException e) {
      // Invalid CIGAR
      Diagnostic.developerLog(e);
      return false;
    }
  }

  @Override
  public int end(final VariantAlignmentRecord rec) {
    return start(rec) + SamUtils.cigarRefLength(rec.getCigar());
  }

  @Override
  public int start(final ReaderRecord<?> var) {
    return var.getStart();
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SamToMatchCigar");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mParams);
    Exam.assertNotNull(mParser);
    Exam.assertTrue(mQdefault >= 0 && mQdefault <= 63);
    return true;
  }
}
