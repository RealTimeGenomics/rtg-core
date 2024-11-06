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
    for (int i = 0; i < mList.size(); ++i) {
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
    for (int i = 0; i < counts; ++i) {
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
