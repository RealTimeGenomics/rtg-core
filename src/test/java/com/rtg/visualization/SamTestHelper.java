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
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, xq);
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
