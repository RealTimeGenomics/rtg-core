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

import java.util.Locale;

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMRecord;

/**
 */
public class SamAssistanceCgLegacy implements SamAssistance {

  private final SamAssistanceSimple mSimple;

  /**
   * Create a new object
   */
  public SamAssistanceCgLegacy() {
    mSimple = new SamAssistanceSimple();
  }

  private static final String INVALID_CIGAR = "Invalid cigar : ";

  static String[] splitCigar(final String cigar, final int split) {
    if (split <= 0) {
      throw new IllegalArgumentException("bad CG split=" + split);
    }
    final StringBuilder sb = new StringBuilder();
    int n = 0;
    int rPos = 0;
    for (int i = 0; i < cigar.length(); ++i) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        assert n > 0;
        switch (c) {
          case SamUtils.CIGAR_SAME_OR_MISMATCH:
          case SamUtils.CIGAR_MISMATCH:
          case SamUtils.CIGAR_SAME:
          case SamUtils.CIGAR_INSERTION_INTO_REF:
          case SamUtils.CIGAR_SOFT_CLIP:
            final int rem = split - rPos;
            assert rem > 0;
            final int todo = n - rem;
            if (todo >= 0) {
              sb.append(rem);
              sb.append(c);
              final String[] res = new String[2];
              res[0] = sb.toString();
              res[1] = (todo > 0 ? "" + todo  + c : "") + cigar.substring(i + 1);
              return res;
            }
            rPos += n;
            break;
          case SamUtils.CIGAR_DELETION_FROM_REF:
          case SamUtils.CIGAR_GAP_IN_READ:
            break;
          default:
            throw new IllegalStateException(INVALID_CIGAR + cigar);
        }
        sb.append(n);
        sb.append(c);
        n = 0;
      }
    }
    throw new IllegalArgumentException("bad CG Cigar=" + cigar + " split=" + split);
  }

  @Override
  public String[] samToReads(final SAMRecord sam, final String template, byte[] templateBytes, final int readStart, final boolean displayDots, boolean displaySoftClip) {
    //System.err.println(template.length());
    if (!sam.hasAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES)) {
      //no overlap use the simple code.
      return mSimple.samToReads(sam, template, templateBytes, readStart, displayDots, displaySoftClip);
    }
    final String gs = sam.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
    final String gc = sam.getStringAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
    final int overLapLength = gs.length() / 2;
    final int overlapPos = Integer.parseInt(gc.substring(0, gc.indexOf('S'))) + overLapLength;
    final String read = sam.getReadString();
    final String cigar = sam.getCigarString();
    final String[] cigars = splitCigar(cigar, overlapPos);
    //System.err.println("gs=" + gs + " gc=" + gc);
    //System.err.println("split cigar " + cigar + " into " + java.util.Arrays.toString(cigars));

    final String[] res = new String[2];
    res[0] = mSimple.cigarToReads(cigars[0], read, template, readStart, displayDots, displaySoftClip);

    final StringBuilder sb = new StringBuilder();
    int overStart = res[0].length();
    for (int i = 0; i < overLapLength; ) {
      if (template.charAt(overStart - 1) != DisplayHelper.INSERT_CHAR) {
        ++i;
      }
      --overStart;
    }
    //System.err.println("readStart: " + readStart + " overstart was: " + length +  " now: " + overStart);
    for (int i = 0, tPos = 0; i < overLapLength; ++tPos) {
      final String te = template.substring(overStart + tPos, overStart + tPos + 1).toUpperCase(Locale.getDefault());
      if (String.valueOf(DisplayHelper.INSERT_CHAR).equals(te)) {
        continue;
      }

      sb.append("1");
      final String ov = String.valueOf(gs.charAt(i + overLapLength));
      if (ov.equals(te)) {
        sb.append(SamUtils.CIGAR_SAME);
      } else {
        sb.append(SamUtils.CIGAR_MISMATCH);
      }
      ++i;
    }
    final String overCigar = sb.toString() + cigars[1];
    final String overRead = gs.substring(overLapLength) + read.substring(overlapPos);
    //System.err.println("overCigar=" + overCigar + " overRead=" + overRead + " overStart=" + overStart);
    res[1] = mSimple.cigarToReads(overCigar, overRead, template, overStart, displayDots, displaySoftClip);
    return res;
  }
}
