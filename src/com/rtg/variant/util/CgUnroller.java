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

import java.util.Arrays;

import com.rtg.reader.CgSamBamSequenceDataSource;
import com.rtg.reader.CgUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.util.Utils;
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
    private final boolean mIsInverted;

    OrientedRead(final byte[] read, final byte[] quality, final boolean isInverted) {
      mRead = read;
      mQuality = quality;
      mIsInverted = isInverted;
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

    public boolean isInverted() {
      return mIsInverted;
    }
  }

  /**
   * Given a SAM record representing a Complete Genomics mapped read, unroll
   * the overlaps and orient the read with the overlap on the left.
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

    final byte[] samRead = Arrays.copyOf(rec.getRead(), rec.getRead().length);
    final byte[] samQualities = rec.getRecalibratedQuality(); // Convert this function to work natively in raw qualities
    final boolean hasQuality = samQualities.length != 0;
    final int samLength = samRead.length;
    if (hasQuality && samLength != samQualities.length) {
      return null;
    }
    final byte[] expandedRead;
    final byte[] expandedQualBytes;
    final String superCigar = rec.getSuperCigar();
    final boolean v1;
    if (superCigar == null) {
      final byte[] gs = rec.getOverlapBases();
      final byte[] gq = rec.getOverlapQuality();
      final String gc = rec.getOverlapInstructions();
      final int overlap = gq.length;
      final boolean legacyLegacy = gq.length == gs.length / 2;
      if (!legacyLegacy && overlap != gs.length) {
        return null;
      }
      expandedRead = CgSamBamSequenceDataSource.unrollLegacyRead(samRead, gs, gc);
      if (expandedRead == null) {
        return null;
      }
      v1 = expandedRead.length == CgUtils.CG_RAW_READ_LENGTH;
      if (!hasQuality) {
        expandedQualBytes = null;
      } else if (legacyLegacy) {
        expandedQualBytes = CgSamBamSequenceDataSource.unrollLegacyLegacyQualities(samQualities, gq, gc);
      } else {
        expandedQualBytes = CgSamBamSequenceDataSource.unrollLegacyRead(samQualities, gq, gc);
      }
    } else {
      final SuperCigarUnroller unroll = new SuperCigarUnroller();
      unroll.setAlignmentRecord(rec);
      unroll.setTemplate(template);
      try {
        unroll.parse();
      } catch (final BadSuperCigarException ex) {
        Diagnostic.developerLog("Ignored SAM CG record due to " + ex.getMessage());
        return null;
      }
      expandedRead = unroll.getByteArray();      //.replaceAll(" ", "");
      v1 = expandedRead.length == CgUtils.CG_RAW_READ_LENGTH;
      final int overlapWidth = unroll.getTemplateOverlapEnd() - unroll.getTemplateOverlapStart();
      final int middle = unroll.getReadOverlapMiddle();
      if (!v1 && expandedRead.length != CgUtils.CG2_RAW_READ_LENGTH && expandedRead.length != CgUtils.CG2_PADDED_LENGTH) {
        return null;
      }
      if (rec.getRead().length != expandedRead.length
        && (overlapWidth <= 0
        //          || samLength + overlapWidth != CG_RAW_READ_LENGTH   NOT necessarily true - deletes in the overlap region change this.
        || middle == -1)) {
        return null;
      }
      if (hasQuality) {
        final int fwdOverlap = v1 ? CgUtils.CG_OVERLAP_POSITION : CgUtils.CG2_OVERLAP_POSITION;
        final int overlapPos = middle <= fwdOverlap ? fwdOverlap : samLength - fwdOverlap;
        final byte[] overlapQual = rec.getOverlapQuality();
        expandedQualBytes = new byte[samQualities.length + overlapQual.length];
        System.arraycopy(samQualities, 0, expandedQualBytes, 0, overlapPos);
        System.arraycopy(overlapQual, 0, expandedQualBytes, overlapPos, overlapQual.length);
        System.arraycopy(samQualities, overlapPos, expandedQualBytes, overlapPos + overlapQual.length, samQualities.length - overlapPos);
//        System.arraycsamQualities.substring(0, overlapPos) + FastaUtils.rawToAsciiString(overlapQual) + samQualities.substring(overlapPos);
        if (expandedQualBytes.length != expandedRead.length) {
          return null;
        }
      } else {
        expandedQualBytes = null;
      }
      //System.out.println("unrollCgRead gives RL=" + rec.getRead().length + " read:" + expandedRead + " (len=" + expandedRead.length() + ")" + "               and qual:" + expandedQual + " (overlapWidth=" + overlapWidth + ")" + " middle=" + middle);
    }
    final boolean cgOverlapOnLeft = v1 && (rec.isFirst() ^ rec.isNegativeStrand()) || !v1 && !rec.isNegativeStrand();
    if (!cgOverlapOnLeft) {
      Utils.reverseInPlace(expandedRead);
      if (expandedQualBytes != null) {
        Utils.reverseInPlace(expandedQualBytes);
      }
    }
    return new CgUnroller.OrientedRead(CgUtils.unPad(expandedRead, !rec.isNegativeStrand()),
      hasQuality ? CgUtils.unPad(expandedQualBytes, !rec.isNegativeStrand()) : null, !cgOverlapOnLeft);
  }

}
