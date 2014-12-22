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
package com.rtg.variant.bayes.complex;

import com.rtg.mode.DnaUtils;
import com.rtg.util.intervals.SequenceNameLocusSimple;

/**
 */
public class ComplexTemplate extends SequenceNameLocusSimple {

  private final byte[] mTemplate;
  private final byte[] mReplaceBytes;
  private final String mReplaceString;

  /**
   * @param template the bytes of the full template.
   * @param refName name of reference sequence.
   * @param start position in the template of the region to be replaced (0 based inclusive).
   * @param end  position in the template of the region to be replaced (0 based exclusive).
   */
  public ComplexTemplate(byte[] template, String refName, int start, int end) {
    super(refName, start, end);
    mTemplate = template;
    final int length = end - start;
    mReplaceBytes = new byte[length];
    final StringBuilder sb = new StringBuilder();
    for (int i = 0, j = start; i < length; i++, j++) {
      final byte b = mTemplate[j];
      mReplaceBytes[i] = b;
      sb.append(DnaUtils.getBase(b));
    }
    mReplaceString = sb.toString();
  }

  /**
   * Get the underlying full template bytes.
   * @return the underlying template bytes.
   */
  public byte[] templateBytes() {
    return mTemplate;
  }

  /**
   * Get the nucleotides in the replacement region as bytes.
   * @return the bytes in the replacement region.
   */
  public byte[] replaceBytes() {
    return mReplaceBytes;
  }

  /**
   * Get the string of nucleotides corresponding to the region to be replaced.
   * @return the string of nucleotides corresponding to the region to be replaced ("" for zero length).
   */
  public String replaceString() {
    return mReplaceString;
  }
}
