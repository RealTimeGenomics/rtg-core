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
    --mReadDeltaPos;
    final byte res;
    if (mIsReverse) {
      final int readPos = mSdfRead.length - mReadPos - 1;
      if (readPos < 0 || readPos >= mSdfRead.length) {
        throw new BadSuperCigarException("Bad reverse simulated super cigar: " + mCigar);
      }
      res = DNA.complement(mSdfRead[readPos]);
      /*
      System.out.println(DnaUtils.bytesToSequenceIncCG(mSdfRead));
      for (int i = 0; i < readPos; ++i) {
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
  protected void doReadOnly(int readNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
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
