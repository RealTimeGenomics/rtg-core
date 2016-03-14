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

package com.rtg.visualization;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMRecord;

/**
 * Given a SAM record, unrolls the reads and stores string representation,
 * represents inserts to template as lowercase character and deletions from template as '-'
 *
 */
class UnrolledRead {

  private final ArrayList<ReadLoc> mList = new ArrayList<>();

  /**
   * Stores an unrolled representation of given read
   *
   * @param rec SAM record to process
   * @param template original template
   * @param zeroBasedTemplateScreenStart  zero based inclusive template start selected on screen (first
   *          location on screen is always 0)
   * @throws BadSuperCigarException on bad cigar
   */
  UnrolledRead(SAMRecord rec, final byte[] template, int zeroBasedTemplateScreenStart) throws BadSuperCigarException {
    //System.err.println("zeroBasedTemplateScreenStart=" + zeroBasedTemplateScreenStart);
    if (isSuperCigar(rec)) {
      final HashMap<Integer, String> unrolledCigar = new SuperCigarUnroller().unroll(rec, template);
      for (Map.Entry<Integer, String> e : unrolledCigar.entrySet()) {
        mList.add(new ReadLoc(e.getKey() + rec.getAlignmentStart() - 1 - zeroBasedTemplateScreenStart, e.getValue()));
      }
    } else if (isLegacyCGCigar(rec)) {
      final HashMap<Integer, String> unrolledCigar = CigarHelper.unrollLegacyCgCigar(rec);
      for (Map.Entry<Integer, String> e : unrolledCigar.entrySet()) {
        mList.add(new ReadLoc(e.getKey() + rec.getAlignmentStart() - 1 - zeroBasedTemplateScreenStart, e.getValue()));
      }
    } else {
      final String unroll = CigarHelper.unrollCigar(rec);
      mList.add(new ReadLoc(rec.getAlignmentStart() - zeroBasedTemplateScreenStart - 1, unroll));
    }
  }

  private boolean isLegacyCGCigar(SAMRecord rec) {
    return rec.hasAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
  }

  private boolean isSuperCigar(SAMRecord rec) {
    return rec.hasAttribute(SamUtils.CG_SUPER_CIGAR);
  }

  /**
   * Returns unrolled reads with correct spacing in screen coordinates
   * @param templateLenghtWithInserts template length including inserts to template
   * @return array of string containing unrolled reads
   */
  String[] lines(int templateLenghtWithInserts) {
    final String[] results = new String[mList.size()];
    for (int i = 0; i < mList.size(); i++) {
      final StringBuilder sb = new StringBuilder();
      final ReadLoc location = mList.get(i);
      //System.err.println("start=" + location.mStart + " end=" + location.mEnd);
      appendSpace(sb, location.mZeroBasedScreenStart);
      sb.append(location.mUnrolledRead);
      appendSpace(sb, templateLenghtWithInserts - sb.length());
      results[i] = sb.toString();
    }
    return results;
  }

  void appendSpace(StringBuilder sb, int counts) {
    for (int i = 0; i < counts; i++) {
      sb.append(' ');
    }
  }

  /**
   * returns if the current location is an insert to the template
   * @param zeroBasedScreenloc 0 based location in screen (console) coordinate, first coordinate is always top-left 0
   * @return true - if its insert, false other wise
   */
  boolean isInsert(int zeroBasedScreenloc) {
    boolean isInsert = false;
    for (ReadLoc rl : mList) {
      if (zeroBasedScreenloc >= rl.mZeroBasedScreenStart && zeroBasedScreenloc < rl.mZeroBasedScreenEnd) {
        final int index = zeroBasedScreenloc - rl.mZeroBasedScreenStart;
        //System.err.println("index=" + index + " mData=" + rl.mData);
        if (Character.isLowerCase(rl.mUnrolledRead.charAt(index))) {
          //System.err.println("XXindex=" + index + " mStart=" + rl.mStart + " mData=" + rl.mData);
          isInsert = true;
        }
      }
    }
    return isInsert;
  }

  /**
   * Holds information regarding unrolled reads
   */
  private static final class ReadLoc {
    int mZeroBasedScreenStart;
    int mZeroBasedScreenEnd;
    String mUnrolledRead;

    ReadLoc(int start, String data) {
      assert start >= 0;
      mZeroBasedScreenStart = start;
      mUnrolledRead = data;
      mZeroBasedScreenEnd = start + data.length();
    }
  }
}
