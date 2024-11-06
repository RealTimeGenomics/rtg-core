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

import java.util.HashMap;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SuperCigarParser;

import htsjdk.samtools.SAMRecord;

/**
 * Handle super cigars in viewer.
 */
public class SuperCigarUnroller extends SuperCigarParser {

  private final StringBuilder mAviewDisplay = new StringBuilder();
  private HashMap<Integer, String> mSnippets = new HashMap<>();
  private int mLastTemplatePosition = 0;
  /** Number of space characters needed for the previous snippet. */
  private Integer mOldSpaces = (Integer) 0;

  private void handleBackwardsMovement() {
    final int p = getTemplatePosition();
    if (p < mLastTemplatePosition) {
      // Going backwards, finish current snippet and start a new one, putting
      // in any necessary spaces to correct for backstep and previous snippet
      mSnippets.put(mOldSpaces, mAviewDisplay.toString());
      final int backShift = mLastTemplatePosition - p + 1;
      mOldSpaces = mAviewDisplay.length() - backShift;
      reset();
    }
    mLastTemplatePosition = p;
  }


  private void reset() {
    mAviewDisplay.setLength(0);
  }


  @Override
  protected void doTemplateSkip(int templateNt) {
    handleBackwardsMovement();
    mAviewDisplay.append(' ');
  }

  @Override
  protected void doReadOnly(int readNt) {
    handleBackwardsMovement();
    mAviewDisplay.append(Character.toLowerCase(DnaUtils.getBase(readNt)));
  }

  @Override
  protected void doTemplateOnly(int templateNt) {
    handleBackwardsMovement();
    mAviewDisplay.append('-');
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    handleBackwardsMovement();
    mAviewDisplay.append(DnaUtils.getBase(readNt));
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    handleBackwardsMovement();
    mAviewDisplay.append(DnaUtils.getBase(nt));
  }

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    doSubstitution(readNt, templateNt);
  }

  @Override
  protected void doUnknownOnRead() {
    handleBackwardsMovement();
    mAviewDisplay.append(DNA.N.name());
  }

  HashMap<Integer, String> unroll(final SAMRecord sam, byte[] templateBytes) throws BadSuperCigarException {
    mSnippets = new HashMap<>();
    reset();
    final String superCigar = sam.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
    // Reset all the parameter information in the parser
    setCigar(superCigar, sam.getStringAttribute(SamUtils.CG_READ_DELTA));
    setTemplate(templateBytes);
    setTemplateStart(sam.getAlignmentStart() - 1);
    mSnippets.clear();
    mLastTemplatePosition = -1;
    // Run the parsing methods generating the AView lines
    parse();
    mSnippets.put(mOldSpaces, mAviewDisplay.toString()); // Handle the last part
    return mSnippets;
  }
}
