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
package com.rtg.variant.util;

import com.rtg.reader.FastaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.VariantAlignmentRecord;

/**
 * Class to unroll CG CIGARs and <code>super CIGARs</code>
 */
public final class CgUnroller {

  private CgUnroller() { }

  /**
   * Hold an oriented read and quality.
   */
  public static class OrientedRead {

    private final byte[] mRead;
    private final byte[] mQuality;
    private final boolean mCgOverlapOnLeft;

    OrientedRead(final byte[] read, final byte[] quality, final boolean cgOverlapOnLeft) {
      mRead = read;
      mQuality = quality;
      mCgOverlapOnLeft = cgOverlapOnLeft;
    }

    public byte[] getRead() {
      return mRead;
    }

    /**
     * Return the quality information or null if no quality available.
     *
     * @return quality
     */
    public byte[] getQuality() {
      return mQuality;
    }

    public boolean isCgOverlapOnLeft() {
      return mCgOverlapOnLeft;
    }
  }

  /**
   * Extract an integer from the start of a string.
   * @param s string containing the integer.
   * @param end last position of integer (0 based exclusive).
   * @return the integer at the start of the string
   */
  private static int stringToInt(final String s, final int start, final int end) {
    int t = 0;
    for (int k = start; k < end; k++) {
      t *= 10;
      t += s.charAt(k) - '0';
    }
    return t;
  }

  private static int nextCigarPos(String cigar, int start) {
    final int end = cigar.length();
    for (int k = start; k < end; k++) {
      switch (cigar.charAt(k)) {
        case 'S':
        case 'G':
          return k;
        default:
      }
    }
    return -1;
  }

  private static String unrollLegacyRead(final String samRead, final String gs, final String gc) {
    final int overlap = gs.length();
    if (overlap == 0) {
      return samRead;
    } else {
      if (gc == null) {
        return null;
      }
      final StringBuilder sbRead = new StringBuilder();
      int lastCigarPos = 0;
      int readPos = 0;
      int attPos = 0;
      while (true) {
        final int cigarPos = nextCigarPos(gc, lastCigarPos);
        if (cigarPos == -1) {
          break;
        }
        final int cigarLen = stringToInt(gc, lastCigarPos, cigarPos);
        if (cigarLen == 0) {
          return null;
        }
        if (gc.charAt(cigarPos) == 'S') {
          sbRead.append(samRead.substring(readPos, readPos + cigarLen));
          lastCigarPos = cigarPos + 1;
          readPos = readPos + cigarLen;
        } else {
          final int consumed = cigarLen * 2;
          if (attPos + consumed > gs.length()) {
            return null;
          }
          sbRead.append(gs.substring(attPos, attPos + consumed));
          attPos = attPos + consumed;
          lastCigarPos = cigarPos + 1;
          readPos += cigarLen;
        }
      }
      if (readPos != samRead.length() || lastCigarPos != gc.length() || attPos != gs.length()) {
        return null;
      }
      return sbRead.toString();
    }
  }

  // This version had a different interpretation of the gq field, where it contained
  // one set of qualities (getting the other from the flattedened representation)
  private static String unrollLegacyLegacyQualities(final String samQualities, final String gq, final String gc) {
    final int overlap = gq.length();
    final int gslength = overlap * 2;
    if (overlap == 0) {
      return samQualities;
    } else {
      if (gc == null) {
        return null;
      }
      final StringBuilder sbQual = new StringBuilder();
      int lastCigarPos = 0;
      int readPos = 0;
      int attPos = 0;
      int att2Pos = 0;
      while (true) {
        final int cigarPos = nextCigarPos(gc, lastCigarPos);
        if (cigarPos == -1) {
          break;
        }
        final int cigarLen = stringToInt(gc, lastCigarPos, cigarPos);
        if (cigarLen == 0) {
          return null;
        }
        if (gc.charAt(cigarPos) == 'S') {
          sbQual.append(samQualities.substring(readPos, readPos + cigarLen));
          lastCigarPos = cigarPos + 1;
          readPos = readPos + cigarLen;
        } else {
          final int consumed = cigarLen * 2;
          if (attPos + consumed > gslength) {
            return null;
          }
          sbQual.append(samQualities.substring(readPos, readPos + cigarLen));
          sbQual.append(gq.substring(att2Pos, att2Pos + cigarLen));
          att2Pos = att2Pos + cigarLen;
          attPos = attPos + consumed;
          lastCigarPos = cigarPos + 1;
          readPos += cigarLen;
        }
      }
      if (readPos != samQualities.length() || lastCigarPos != gc.length() || attPos != gslength) {
        return null;
      }
      return sbQual.toString();
    }
  }

  /**
   * Given a SAM record representing a Complete Genomics mapped read, unroll
   * the overlaps are orient the read with the large gap on the right.
   * If the input record is found to be invalid null is returned.
   *
   * @param rec SAM record
   * @param template template bytes
   * @return unrolled representation
   */
  public static CgUnroller.OrientedRead unrollCgRead(final VariantAlignmentRecord rec, final byte[] template) {
    // Produce a read exactly 35 nucleotides long, making it look
    // like a forward-complement left-arm.
    if (!rec.isReadPaired()) {
      // Can't determine arm without this information
      return null;
    }

    final String samRead = new String(rec.getRead());
    final String samQualities = new String(FastaUtils.rawToAsciiQuality(rec.getQuality())); // Convert this function to work natively in raw qualities
    final boolean hasQuality = samQualities.length() != 0;
    final int samLength = samRead.length();
    if (hasQuality && samLength != samQualities.length()) {
      return null;
    }
    final String expandedRead;
    final String expandedQual;
    final String superCigar = rec.getSuperCigar();
    if (superCigar == null) {
      final String gs = rec.getOverlapBases();
      final String gq = rec.getOverlapQuality();
      final String gc = rec.getOverlapInstructions();
      final int overlap = gq.length();
      final boolean legacyLegacy = gq.length() == gs.length() / 2;
      if (!legacyLegacy && overlap != gs.length()) {
        return null;
      }
      expandedRead = unrollLegacyRead(samRead, gs, gc);
      if (expandedRead == null) {
        return null;
      }
      if (!hasQuality) {
        expandedQual = null;
      } else if (legacyLegacy) {
        expandedQual = unrollLegacyLegacyQualities(samQualities, gq, gc);
      } else {
        expandedQual = unrollLegacyRead(samQualities, gq, gc);
      }
    } else {
      final SuperCigarUnroller unroll = new SuperCigarUnroller();
      unroll.setAlignmentRecord(rec);
      unroll.setTemplate(template, template.length);
      try {
        unroll.parse();
      } catch (final BadSuperCigarException ex) {
        Diagnostic.developerLog("Ignored SAM CG record due to " + ex.getMessage());
        return null;
      }
      expandedRead = unroll.getString();      //.replaceAll(" ", "");
      final int overlapWidth = unroll.getTemplateOverlapEnd() - unroll.getTemplateOverlapStart();
      final int middle = unroll.getReadOverlapMiddle();
      if (hasQuality) {
        final int overlapPos = middle <= SamUtils.CG_OVERLAP_POSITION ? SamUtils.CG_OVERLAP_POSITION : samLength - SamUtils.CG_OVERLAP_POSITION;
        final String overlapQual = rec.getOverlapQuality();
        expandedQual = samQualities.substring(0, overlapPos) + overlapQual + samQualities.substring(overlapPos);
      } else {
        expandedQual = null;
      }
      //System.out.println("unrollCgRead gives RL=" + rec.getRead().length + " read:" + expandedRead + " (len=" + expandedRead.length() + ")" + "               and qual:" + expandedQual + " (overlapWidth=" + overlapWidth + ")" + " middle=" + middle);
      if (rec.getRead().length != SamUtils.CG_RAW_READ_LENGTH && //if the read in the sam record is the normal cg read length (i.e. no overlap), these further checks are invalid
          (overlapWidth <= 0
           || expandedRead.length() != SamUtils.CG_RAW_READ_LENGTH
           //          || samLength + overlapWidth != CG_RAW_READ_LENGTH   NOT necessarily true - deletes in the overlap region change this.
           || middle == -1
           || (expandedQual != null && expandedQual.length() != SamUtils.CG_RAW_READ_LENGTH))) {
        return null;
      }
    }
    final boolean cgOverlapOnLeft = rec.isCgOverlapLeft();
    final String reversedRead;
    final String reversedQual;
    if (cgOverlapOnLeft) {
      reversedRead = expandedRead;
      reversedQual = expandedQual;
    } else {
      reversedRead = StringUtils.reverse(expandedRead);
      reversedQual = StringUtils.reverse(expandedQual);
    }
    return new CgUnroller.OrientedRead(reversedRead.getBytes(), hasQuality ? FastaUtils.asciiToRawQuality(reversedQual) : null, cgOverlapOnLeft);
  }

}
