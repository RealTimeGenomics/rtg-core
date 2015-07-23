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

import com.rtg.alignment.CgGotohEditDistance;
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
      assert read.length == CgGotohEditDistance.CG_RAW_READ_LENGTH;
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
  private static int stringToInt(final String s, final int end) {
    int t = 0;
    for (int k = 0; k < end; k++) {
      t *= 10;
      t += s.charAt(k) - '0';
    }
    return t;
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
      final int overlap = gq.length();
      if (2 * overlap != gs.length()) {
        return null;
      }
      if (overlap + samLength != CgGotohEditDistance.CG_RAW_READ_LENGTH) {
        return null;
      }
      if (overlap == 0) {
        expandedRead = samRead;
        expandedQual = hasQuality ? samQualities : null;
      } else {
        final int overlapPos;
        final String gc = rec.getOverlapInstructions();
        if (gc == null) {
          return null;
        }
        overlapPos = stringToInt(gc, gc.indexOf('S'));
        if (overlapPos == 0) {
          return null;
        }
        expandedRead = samRead.substring(0, overlapPos) + gs + samRead.substring(overlapPos + overlap);
        expandedQual = !hasQuality ? null : samQualities.substring(0, overlapPos + overlap) + gq + samQualities.substring(overlapPos + overlap);
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
