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
package com.rtg.sam;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;


/**
 * Validates a simulated super cigar against the template and original read
 */
public class SimulatedSuperCigarUnroller extends SuperCigarParser {
  private byte[] mSdfRead = null;
  protected StringBuilder mBuilder;
  private boolean mIsReverse = false;
  //private byte[] mFlattenedRead = null; //pos specified by mReadPos in super
  //private int mFlattenedReadPos;

  //private int mMismatches;

  public String getString() {
    return mBuilder.toString();
  }

  @Override
  public void setCigar(String superCigar, String readDelta) {
    if (superCigar.contains(":")) {
      throw new RuntimeException("old style simulated cigar found: " + superCigar);
    }
    super.setCigar(superCigar, readDelta);
  }

  /**
   * Set read as it appears in original SDF
   * @param sdfRead the read data
   */
  public void setSdfRead(byte[] sdfRead) {
    mSdfRead = sdfRead;
  }

  /**
   * Set read as it appears in original SDF
   * @param sdfRead the read data
   */
  public void setSdfRead(String sdfRead) {
    setSdfRead(DnaUtils.encodeString(sdfRead));
  }

  public void setReverse(boolean reverse) {
    mIsReverse = reverse;
  }

  private void reset() {
    mBuilder = new StringBuilder();
  }

  /**
   * @throws BadSuperCigarException on bad cigar.
   */
  @Override
  public void parse() throws BadSuperCigarException {
    reset();
    super.parse();
  }

  @Override
  protected byte getReadDelta(int pos) throws BadSuperCigarException {
    // We don't actually have a read delta so improvise
    mReadDeltaPos--;
    final byte res;
    if (mIsReverse) {
      final int readPos = mSdfRead.length - mReadPos - 1;
      if (readPos < 0 || readPos >= mSdfRead.length) {
        throw new BadSuperCigarException("Bad reverse simulated super cigar: " + mCigar);
      }
      res = DNA.complement(mSdfRead[readPos]);
      /*
      System.out.println(DnaUtils.bytesToSequenceIncCG(mSdfRead));
      for (int i = 0; i < readPos; i++) {
        System.out.print(" ");
      }
        System.out.println("^");
       *
       */
    } else {
      if (mReadPos < 0 || mReadPos >= mSdfRead.length) {
        throw new BadSuperCigarException("Bad simulated super cigar: " + mCigar);
      }
       res = mSdfRead[mReadPos];
    }
    return res;
  }


  @Override
  protected void doReadHardClip() {
  }

  @Override
  protected void doReadSoftClip(int readNt) {
  }

  @Override
  protected void doTemplateOverlap() {
  }

  @Override
  protected void doTemplateSkip(int templateNt) {
  }

  @Override
  protected void doReadOnly(int readNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doTemplateOnly(int templateNt) {
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    mBuilder.append(DNA.valueChars()[nt]);
  }

  @Override
  protected void advanceTemplate(boolean direction) {
      super.advanceTemplate(direction);
  }

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doUnknownOnRead() {
    mBuilder.append(DNA.N.name());
  }

}
