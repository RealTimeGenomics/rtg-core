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

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 */
final class SamTestHelper {



  private SamTestHelper() { }

  static final String TEMPLATE_NAME = "simulatedSequence1";
  static final String READ_ID = "0";

  static SAMRecord getSuperCigarSAMRecord(String cigar, String read, String xs, String xq, String xr) {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setAlignmentStart(1);
    rec.setCigarString(cigar);
    rec.setReadName(READ_ID);
    rec.setReferenceName(TEMPLATE_NAME);
    rec.setReadBases(read.getBytes());
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, xs);
    rec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, xq);
    rec.setAttribute(SamUtils.CG_READ_DELTA, xr);
    return rec;
  }

  static SAMRecord getLegacyCGSamRecord(String read, String cigar, String gs, String gc, String gq) {
    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    r.setAlignmentStart(1);
    r.setCigarString(cigar);
    r.setReadName(READ_ID);
    r.setReferenceName(TEMPLATE_NAME);
    r.setReadBases(read.getBytes());
    r.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, gs);
    r.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, gc);
    r.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, gq);
    return r;
  }

  static SAMRecord getSamRecord(String read, String cigar) {
    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    r.setAlignmentStart(1);
    r.setCigarString(cigar);
    r.setReadName(READ_ID);
    r.setReferenceName(TEMPLATE_NAME);
    r.setReadBases(read.getBytes());
    return r;
  }
}
