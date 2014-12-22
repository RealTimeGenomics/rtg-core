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
package com.rtg.util.intervals;

import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.Locale;

/**
 * Class to handle single-sequence region restrictions on genome position data (SNP, BED, SAM)
 * String representations use values as 1-based with the end inclusive. Allows start
 * and end values to be undefined to indicate regions that extend to the end or
 * encompass an entire sequence. Internal API calls are all zero-based half-open.
 *
 */
public class RegionRestriction implements SequenceNameLocus {

  /** Used to indicate start or end is missing */
  public static final int MISSING = -1;

  private final String mSequence;
  private final int mStart; // 0-based inclusive
  private final int mEnd;   // 0-based exclusive (same as 1-based inclusive)

  /**
   * Restriction to named sequence, accepts restrictions of
   * the form name, name:start, name:start-end or name:start+length
   * @param restriction restriction string
   */
  public RegionRestriction(String restriction) {
    if (restriction.length() == 0) {
      throw malformedRange(restriction);
    }
    final int rangepos = restriction.lastIndexOf(':');
    if (rangepos != -1) {
      mSequence = restriction.substring(0, rangepos);
      final String range = restriction.substring(rangepos + 1);
      final NumberFormat fmt = NumberFormat.getInstance(Locale.getDefault());
      fmt.setParseIntegerOnly(true);
      final boolean useLength = range.contains("+");
      final String[] parts;
      if (useLength) {
        if (range.contains("-")) {
          throw malformedRange(restriction);
        }
        parts = range.split("\\+");
      } else {
        parts = range.split("-");
      }
      if (parts.length != 2 && parts.length != 1) {
        throw malformedRange(restriction);
      }
      final ParsePosition pos = new ParsePosition(0);
      final Number start = fmt.parse(parts[0], pos);
      if (parts[0].length() != pos.getIndex() || pos.getIndex() == 0) {
        throw malformedRange(restriction);
      }
      mStart = start.intValue() - 1;
      if (parts.length == 1) {
        mEnd = MISSING;
      } else {
        pos.setIndex(0);
        final Number last = fmt.parse(parts[1], pos);
        if (parts[1].length() != pos.getIndex() || pos.getIndex() == 0) {
          throw malformedRange(restriction);
        }
        mEnd = (useLength ? mStart : 0) + last.intValue();
      }
      if ((mStart < 0) || (mEnd != MISSING && mEnd < (mStart + 1))) {
        throw malformedRange(restriction);
      }
    } else {
      mSequence = restriction;
      mStart = MISSING;
      mEnd = MISSING;
    }
  }

  private IllegalArgumentException malformedRange(String restriction) {
    return new IllegalArgumentException("Malformed range in restriction: \"" + restriction + "\"");
  }

  /**
   * Restriction to region in named sequence
   * @param sequence sequence name
   * @param start start position of allowed region (0 based inclusive)
   * @param end end position of allowed region (0 based exclusive)
   */
  public RegionRestriction(String sequence, int start, int end) {
    if (sequence == null) {
      throw new NullPointerException();
    }
    mSequence = sequence;
    mStart = start;
    mEnd = end;
  }

  /**
   * @return sequence name for restriction
   */
  @Override
  public String getSequenceName() {
    return mSequence;
  }

  /**
   * @return start position for restriction (0 based inclusive).
   */
  @Override
  public int getStart() {
    return mStart;
  }

  /**
   * @return end position for restriction (0 based exclusive)
   */
  @Override
  public int getEnd() {
    return mEnd;
  }

  @Override
  public int getLength() {
    return (mStart == MISSING || mEnd == MISSING) ? -1 : mEnd - mStart;
  }

  /**
   * @return string of form <code>sequenceName[:startPos[-endPos]]</code>
   */
  @Override
  public String toString() {
    return mSequence + (mStart != MISSING ? ":" + (mStart + 1) + (mEnd != MISSING ? "-" + mEnd : "") : "");
  }

  /**
   * Method to validate the region provided.
   * @param region the region string to test
   * @return true if the string conforms to the parsing rules, false otherwise
   */
  public static boolean validateRegion(String region) {
    try {
      new RegionRestriction(region);
    } catch (final IllegalArgumentException e) {
      return false;
    }
    return true;
  }

  /**
   * Test whether the specified range overlaps the current region
   * @param locus the sequence of interest
   * @return true if the region overlaps the supplied parameters
   */
  @Override
  public boolean overlaps(SequenceNameLocus locus) {
    return overlaps(locus.getSequenceName(), locus.getStart(), locus.getEnd());
  }

  /**
   * Test whether the specified range overlaps the current region
   * @param sequence the sequence of interest
   * @param start the start of the region (0-based, inclusive)
   * @param end the end of the region (0-based, exclusive)
   * @return true if the region overlaps the supplied parameters
   */
  public boolean overlaps(String sequence, int start, int end) {
    if (start == MISSING || end == MISSING || end < start) {
      throw new IllegalArgumentException();
    }
    if (!mSequence.equals(sequence)) {
      return false;
    }
    if ((mEnd != MISSING && start >= mEnd)
        || (end <= mStart)) {
      return false;
    }
    return true;
  }

  @Override
  public boolean contains(String sequence, int pos) {
    return SequenceNameLocusSimple.contains(this, sequence, pos);
  }

}
