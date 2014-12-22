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

import java.util.Arrays;

/**
 * A mock SAM/BAM record.
 * Usage: call setField or addAttribute for each field your test needs to read.
 *   You can also reuse a record, by overwriting just the fields that change.
 *
 */
public class MockSamBamRecord implements SamBamRecord {

  // stores all fields, by their position.
  private String[] mField = new String[SamBamConstants.ATTRIBUTES_FIELD];

  /**
   * @param attr a string like "NM:i:32".
   */
  public void addAttribute(String attr) {
    assert attr.matches("[A-Z][A-Z]:[AifZH]:.+");
    mField = Arrays.copyOf(mField, mField.length + 1);
    mField[mField.length - 1] = attr;
  }

  @Override
  public String getReadName() {
    return mField[SamBamConstants.RNAME_FIELD];
  }

  @Override
  public int getReadNumber() {
    return Integer.parseInt(mField[SamBamConstants.RNAME_FIELD]);
  }

  @Override
  public int getFlags() {
    return Integer.parseInt(mField[SamBamConstants.FLAG_FIELD]);
  }

  @Override
  public boolean hasAttribute(String tag) {
    return getFieldNumFromTag(tag) >= 0;
  }

  @Override
  public char getAttributeType(String tag) {
    final int pos = getFieldNumFromTag(tag);
    return pos < 0 ? 0 : mField[pos].charAt(3);
  }

  @Override
  public int getFieldNumFromTag(String tag) {
    for (int i = SamBamConstants.ATTRIBUTES_FIELD; i < mField.length; i++) {
      if (mField[i] != null && mField[i].startsWith(tag)) {
        return i;
      }
    }
    return -1;
  }

  @Override
  public String[] getAttributeTags() {
    throw new UnsupportedOperationException();
  }

  @Override
  public Object getAttributeValue(String tag) {
    final int pos = getFieldNumFromTag(tag);
    return pos < 0 ? null : mField[pos].substring(5);
  }

  @Override
  public int getIntAttribute(String tag) {
    final int pos = getFieldNumFromTag(tag);
    return pos < 0 ? 0 : Integer.parseInt(mField[pos].substring(5));
  }

  @Override
  public int getNumFields() {
    return mField.length;
  }

  @Override
  public String getField(int fieldNum) {
    return mField[fieldNum];
  }

  @Override
  public int getIntField(int fieldNum) {
    return Integer.parseInt(mField[fieldNum]);
  }

  @Override
  public byte[] getFieldBytes(int fieldNum) {
    return mField[fieldNum].getBytes();
  }

}
