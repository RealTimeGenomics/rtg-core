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
    for (int i = 0; i < cigar.length(); i++) {
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
  public String[] samToReads(final SAMRecord sam, final String template, byte[] templateBytes, final int readStart, final boolean displayDots) {
    //System.err.println(template.length());
    if (!sam.hasAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES)) {
      //no overlap use the simple code.
      return mSimple.samToReads(sam, template, templateBytes, readStart, displayDots);
    }
    final String gs = sam.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);
    assert gs.length() > 0 && (gs.length() & 1) == 0;
    final String gc = sam.getStringAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
    final String read = sam.getReadString();
    final String cigar = sam.getCigarString();
    //System.err.println("sam=" + sam.getReadName() + " " + sam.format());
    final int overLapLength = gs.length() >> 1;
    final int overlapPos = Integer.parseInt(gc.substring(0, gc.indexOf('S'))) + overLapLength;
    final String[] cigars = splitCigar(cigar, overlapPos);
    //System.err.println("gs=" + gs + " gc=" + gc);
    //System.err.println("split cigar " + cigar + " into " + java.util.Arrays.toString(cigars));

    final String[] res = new String[2];
    res[0] = mSimple.cigarToReads(cigars[0], read, template, readStart, displayDots);

    final StringBuilder sb = new StringBuilder();
    int overStart = res[0].length();
    for (int i = 0; i < overLapLength; ) {
      if (template.charAt(overStart - 1) != '_') {
        i++;
      }
      overStart--;
    }
    //System.err.println("readStart: " + readStart + " overstart was: " + length +  " now: " + overStart);
    for (int i = 0, tPos = 0; i < overLapLength; tPos++) {
      final String te = template.substring(overStart + tPos, overStart + tPos + 1).toUpperCase(Locale.getDefault());
      if (te.equals("_")) {
        continue;
      }

      sb.append("1");
      final String ov = "" + gs.charAt(i + overLapLength);
      if (ov.equals(te)) {
        sb.append(SamUtils.CIGAR_SAME);
      } else {
        sb.append(SamUtils.CIGAR_MISMATCH);
      }
      i++;
    }
    final String overCigar = sb.toString() + cigars[1];
    final String overRead = gs.substring(overLapLength) + read.substring(overlapPos);
    //System.err.println("overCigar=" + overCigar + " overRead=" + overRead + " overStart=" + overStart);
    res[1] = mSimple.cigarToReads(overCigar, overRead, template, overStart, displayDots);
    return res;
  }
}
