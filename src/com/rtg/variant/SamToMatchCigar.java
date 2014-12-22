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

package com.rtg.variant;

import java.util.TreeSet;

import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.CigarFormatter;
import com.rtg.sam.ReaderRecord;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
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
  public boolean process(final byte[] templateBytes, final VariantAlignmentRecord rec, final TreeSet<VariantAlignmentRecord> samComplexContext) {
    final AbstractMachineErrorParams me = mChooser.machineErrors(rec);
    final int readScore = VariantUtils.readScoreFromAlignmentRecord(rec, mParams);
    if (readScore < 0) {
      // Invalid read score
      return false;
    }
    try {
      // System.err.println("read=" + read + " cigar=" + cigarString);
      mParser.toMatcher(me, rec, mQdefault, templateBytes);
      if (samComplexContext != null) {
        samComplexContext.add(rec);
      }
      return true;
    } catch (final BadSuperCigarException e) {
      // Invalid CIGAR
      Diagnostic.developerLog(e);
      return false;
    }
  }

  @Override
  public int end(final VariantAlignmentRecord rec) {
    return start(rec) + CigarFormatter.cigarRefLength(rec.getCigar());
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
