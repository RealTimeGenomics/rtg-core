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
package com.rtg.simulation;

import com.rtg.util.StringUtils;

/**
 * Parses names generated by the old read simulator.
 */
public class OldReadNameParser implements SimulatedReadNameParser {

  static final char READ_SIM_SEPARATOR = ':';
  static final int EXPECT_MIN_FIELDS = 6;

  /** Index of the read label field. */
  public static final int READ_LABEL = 0;
  /** Index of the template sequence name field. */
  public static final int TEMPLATE_SEQ_NAME = 1;
  /** Index of the template position/frame combo field. */
  public static final int TEMPLATE_POS_ID = 2;
  /** Index of the substitution count field. */
  public static final int SUBS_COUNT_ID = 3;
  /** Index of the insertion count field. */
  public static final int INS_COUNT_ID = 4;
  /** Index of the deletion count field. */
  public static final int DEL_COUNT_ID = 5;

  /**
   * @param readName an example read name
   * @return true if this read name looks like something this parser deals with
   */
  public static boolean looksOk(String readName) {
    if (readName.indexOf(READ_SIM_SEPARATOR) != -1) {
      final String[] f = StringUtils.split(readName, READ_SIM_SEPARATOR);
      if (f.length >= EXPECT_MIN_FIELDS) {
        return true;
      }
    }
    return false;
  }



  private String[] mFields = null;
  private int mId = 0;
  private int mReadLength = 0;
  private boolean mReverse = false;
  private int mSubs;
  private int mIns;
  private int mDel;


  @Override
  public boolean setReadInfo(String readName, int readLength) {
    mFields = StringUtils.split(readName, READ_SIM_SEPARATOR);
    if (mFields.length < 6) {
      return false;
    }
    mId = Integer.parseInt(mFields[READ_LABEL].replaceAll("read", ""));
    mReadLength = readLength;
    final String pos = mFields[TEMPLATE_POS_ID];
    mReverse = pos.charAt(pos.length() - 1) == 'R';
    if (mReverse) {
      mFields[TEMPLATE_POS_ID] = pos.substring(0, pos.length() - 1);
    }
    mSubs = Integer.parseInt(mFields[SUBS_COUNT_ID].replaceAll("S", ""));
    mIns = Integer.parseInt(mFields[INS_COUNT_ID].replaceAll("I", ""));
    mDel = Integer.parseInt(mFields[DEL_COUNT_ID].replaceAll("D", ""));
    return true;
  }

  @Override
  public boolean isChimera() {
    return false;
  }

  @Override
  public boolean isDuplicate() {
    return false;
  }

  @Override
  public long templatePosition() {
    final long pos = Long.parseLong(mFields[TEMPLATE_POS_ID]);
    return mReverse ? pos - mReadLength : pos;
  }

  @Override
  public String templateName() {
    return mFields[TEMPLATE_SEQ_NAME];
  }

  @Override
  public int templateSet() {
    return 0;
  }

  @Override
  public boolean forwardFrame() {
    return !mReverse;
  }

  @Override
  public int readId() {
    return mId;
  }

  @Override
  public String readName() {
    return mFields[READ_LABEL];
  }

  @Override
  public int substitutions() {
    return mSubs;
  }

  @Override
  public int insertions() {
    return mIns;
  }

  @Override
  public int deletions() {
    return mDel;
  }

  @Override
  public int numMismatches() {
    return mSubs + mIns + mDel;
  }
}
