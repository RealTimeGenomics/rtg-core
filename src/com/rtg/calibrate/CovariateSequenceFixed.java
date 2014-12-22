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

package com.rtg.calibrate;

import java.util.List;
import java.util.Locale;

import net.sf.samtools.SAMRecord;

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
