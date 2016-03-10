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

import java.util.ArrayList;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SuperCigarParser;

import htsjdk.samtools.SAMRecord;

/**
 * Handle super cigars in viewer.
 */
public class SamAssistanceCg extends SuperCigarParser implements SamAssistance {

  private final StringBuilder mAviewDisplay = new StringBuilder();
  private final ArrayList<String> mSnippets = new ArrayList<>();
  private String mExpandedTemplate = null;
  private int mExpandedTemplatePosition = 0;
  private int mLastTemplatePosition = 0;
  private boolean mWasInsert = false;
  private boolean mDisplayDots = true; //on by default
  /** Number of space characters needed for the previous snippet. */
  private int mOldSpaces = 0;

  private void resetSnippet() {
    mAviewDisplay.setLength(0);
    for (int k = 0; k < mOldSpaces; k++) {
      mAviewDisplay.append(' ');
    }
  }

  private void handleBackwardsMovement() {
    final int p = getTemplatePosition();
    if (p <= mLastTemplatePosition && !mWasInsert) {
      // Going backwards, finish current snippet and start a new one, putting
      // in any necessary spaces to correct for backstep and previous snippet
      mSnippets.add(mAviewDisplay.toString());
      final int backShift = mLastTemplatePosition - p + 1;
      int extras = 0;
      for (int k = 0; k < backShift; k++) {
        while (--mExpandedTemplatePosition >= 0
               && mExpandedTemplate.charAt(mExpandedTemplatePosition) == '_') {
          extras++;
        }
      }
      mOldSpaces = mAviewDisplay.length() - backShift - extras;
      resetSnippet();
    }
    mWasInsert = false;
    mLastTemplatePosition = p;
  }

  private void handleExpandedTemplate(final int shift) {
    int s = shift;
    while (mExpandedTemplatePosition < mExpandedTemplate.length()
           && mExpandedTemplate.charAt(mExpandedTemplatePosition) == '_') {
      if (s == 0) {
        mAviewDisplay.append('_');
      } else {
        s--;
      }
      mExpandedTemplatePosition++;
    }
    mExpandedTemplatePosition += 1 - shift;
  }
  private void consumeInsert() {
    mExpandedTemplatePosition++;
  }

  @Override
  protected void doTemplateSkip(int templateNt) {
    handleBackwardsMovement();
    handleExpandedTemplate(0);
    mAviewDisplay.append(' ');
  }

  @Override
  protected void doReadOnly(int readNt) {
    handleBackwardsMovement();
    mAviewDisplay.append(DnaUtils.getBase(readNt));
    consumeInsert();
    mWasInsert = true;
  }

  @Override
  protected void doTemplateOnly(int templateNt) {
    handleBackwardsMovement();
    handleExpandedTemplate(0);
    mAviewDisplay.append('-');
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    handleBackwardsMovement();
    handleExpandedTemplate(0);
    mAviewDisplay.append(DnaUtils.getBase(readNt));
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    handleBackwardsMovement();
    handleExpandedTemplate(0);
    if (mDisplayDots) {
      mAviewDisplay.append('.');
    } else {
      mAviewDisplay.append(DnaUtils.getBase(nt));
    }
  }

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    doSubstitution(readNt, templateNt);
  }

  @Override
  protected void doUnknownOnRead() {
    handleBackwardsMovement();
    handleExpandedTemplate(0);
    mAviewDisplay.append('N');
  }

  @Override
  public String[] samToReads(final SAMRecord sam, final String template, byte[] templateBytes, final int readStart, final boolean displayDots) throws BadSuperCigarException {
    final String superCigar = sam.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
    // Reset all the parameter information in the parser
    setCigar(superCigar, sam.getStringAttribute(SamUtils.CG_READ_DELTA));
    setTemplate(templateBytes);
    setTemplateStart(readStart);
    mSnippets.clear();
    mLastTemplatePosition = -1;
    mExpandedTemplate = template;
    mExpandedTemplatePosition = readStart;
    mOldSpaces = readStart;
    mDisplayDots = displayDots;
    resetSnippet();

    // Run the parsing methods generating the AView lines
    parse();
    mSnippets.add(mAviewDisplay.toString()); // Handle the last part
    return mSnippets.toArray(new String[mSnippets.size()]);
  }
}
