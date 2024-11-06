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
package com.rtg.variant.util;

import java.io.ByteArrayOutputStream;

import com.rtg.mode.DNA;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SuperCigarParser;
import com.rtg.variant.VariantAlignmentRecord;

/**
 * Validates a super cigar against the template and original read.
 * This reconstructs the original read from the alignment record plus the template.
 * Usage example:
 * <pre>
 *   SuperCigarUnroller unroll = new SuperCigarUnroller();
 *   unroll.setSamRecord(samrec);
 *   unroll.setTemplate(template);
 *   unroll.parse();
 *   String originalRead = unroll.getString();
 * </pre>
 *
 */
public class SuperCigarUnroller extends SuperCigarParser {

  protected ByteArrayOutputStream mBaos = new ByteArrayOutputStream();

  /**
   * Set the SAM record, which must contain a super cigar and read delta attribute.
   * @param alignmentRecord a valid SAM record.
   */
  public void setAlignmentRecord(VariantAlignmentRecord alignmentRecord) {
    final String superCigar = alignmentRecord.getSuperCigar();
    final String cigar =  superCigar == null ? alignmentRecord.getCigar() : superCigar;
    super.setCigar(cigar, alignmentRecord.getCGReadDelta());
    setTemplateStart(alignmentRecord.getStart());
  }


  /**
   * @return the read string, with spaces for CG gaps.
   */
  public byte[] getByteArray() {
    return mBaos.toByteArray();
  }

  /**
   * @throws BadSuperCigarException on bad cigar.
   */
  @Override
  public void parse() throws BadSuperCigarException {
    mBaos = new ByteArrayOutputStream();
    super.parse();
    //checkMismatchCounts();
  }

  @Override
  protected void doReadSoftClip(int readNt) {
    mBaos.write(readNt);
  }

  @Override
  protected void doReadOnly(int readNt) {
    mBaos.write(readNt);
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    mBaos.write(readNt);
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    mBaos.write(nt);
  }

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    mBaos.write(readNt);
  }

  @Override
  protected void doUnknownOnRead() {
    mBaos.write(DNA.N.ordinal());
  }
}
