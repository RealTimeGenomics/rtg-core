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

package com.rtg.bed;

import static com.rtg.util.StringUtils.TAB;

import java.util.Arrays;

import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.util.StringUtils;

/**
 * A BED record class.  Positions in BED format are zero-based.
 *
 */
public class BedRecord extends SequenceNameLocusSimple {

  protected final String[] mAnnotations;

  /**
   * Constructor.
   * @param sequence the sequence name
   * @param start the start of the region (0-based inclusive)
   * @param end the end of the region (0-based exclusive)
   * @param annotations the list of annotation fields
   */
  public BedRecord(String sequence, int start, int end, String... annotations) {
    super(sequence, start, end);
    if (annotations == null) {
      mAnnotations = new String[0];
    } else {
      mAnnotations = annotations;
    }
  }

  /**
   * Get annotations, which are from columns 3 onwards in the BED file.
   * @return Returns the annotations.
   */
  public String[] getAnnotations() {
    return mAnnotations;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(mSequence);
    sb.append(TAB).append(getStart());
    sb.append(TAB).append(getEnd());
    for (final String s : mAnnotations) {
      sb.append(TAB).append(s);
    }
    return sb.toString();
  }

  /**
   * Get a BedRecord from a single line.
   * @param line the BED line to parse
   * @return the BedRecord
   * @throws NumberFormatException if the start or end cannot be parsed
   * @throws ArrayIndexOutOfBoundsException if there are not enough fields
   */
  public static BedRecord fromString(String line) throws NumberFormatException, ArrayIndexOutOfBoundsException {
    final String[] parts = StringUtils.split(line, '\t');
    return new BedRecord(parts[0], Integer.parseInt(parts[1]), Integer.parseInt(parts[2]), Arrays.copyOfRange(parts, 3, parts.length));
  }
}
