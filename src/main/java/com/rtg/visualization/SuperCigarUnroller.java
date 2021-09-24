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
