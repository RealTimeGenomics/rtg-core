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
package com.rtg.tabix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Used as basis for reading files for <code>TABIX</code> indexing
 */
@TestClass("com.rtg.tabix.SamPositionReaderTest")
abstract class AbstractPositionReader implements BlockCompressedPositionReader {

  private final int[] mTabs;
  private int mTabsUsed;
  private String mCurrentLine;
  private final BlockCompressedLineReader mReader;
  private final char mMeta;
  private String mNextLine;
  private long mNextVirtualOffset;
  private long mCurrentVirtualOffset;

  protected int mRefId;
  protected String mLastReferenceName;
  protected String mReferenceName;
  protected int mStartPosition;
  protected int mLengthOnReference;
  protected int mBin;
  protected HashMap<String, Integer> mSequenceNames;
  protected ArrayList<String> mNamesList;

  /**
   * Constructor
   * @param reader source of file
   * @param numColumns number of columns in file
   * @param meta character signifying start of header line
   * @param skip number of lines to skip at beginning of file
   */
  AbstractPositionReader(BlockCompressedLineReader reader, int numColumns, char meta, int skip) {
    mReader = reader;
    mTabs = new int[numColumns + 1];
    mTabs[0] = -1;
    mTabsUsed = 1;
    mRefId = -1;
    mMeta = meta;
    mSequenceNames = new HashMap<>();
    mNamesList = new ArrayList<>();
    for (int i = 0; i < skip; i++) {
      mReader.readLine();
    }
    populateNext();
  }

  @Override
  public String getRecord() {
    return mCurrentLine;
  }

  @Override
  public void seek(long virtualOffset) throws IOException {
    mReader.seek(virtualOffset);
    mCurrentLine = null;
    mCurrentVirtualOffset = 0;
    mNextVirtualOffset = 0;
    mStartPosition = 0;
    mLengthOnReference = 0;
    mBin = 0;
    mReferenceName = null;
    mLastReferenceName = null;
    mRefId = Integer.MIN_VALUE;
    mTabsUsed = 1;
    populateNext();
  }

  /**
   * implementors should set the {@link AbstractPositionReader#mReferenceName} field to
   * the name of the reference sequence that applies to the current record
   */
  protected abstract void setReferenceName() throws IOException;

  /**
   * implementors should set the {@link AbstractPositionReader#mStartPosition} field to
   * the start position (0-based) on the reference of the current record and the {@link AbstractPositionReader#mLengthOnReference}
   * field to the length of the region the current record applies to
   */
  protected abstract void setStartAndLength() throws IOException;

  @Override
  public void next() throws IOException {

    mCurrentLine = mNextLine;
    mCurrentVirtualOffset = mNextVirtualOffset;
    mTabsUsed = 1;
    populateNext();
    mLastReferenceName = mReferenceName;
    setReferenceName();
    setStartAndLength();
    mBin = TabixIndexer.reg2bin(mStartPosition, mStartPosition + mLengthOnReference);
    if (mLastReferenceName == null || !mLastReferenceName.equals(mReferenceName)) {
      final Integer pos = mSequenceNames.get(mReferenceName);
      if (pos != null) {
        mRefId = pos;
      } else {
        mRefId = mSequenceNames.size();
        mSequenceNames.put(mReferenceName, mSequenceNames.size());
        mNamesList.add(mReferenceName);
      }
    }
  }

  @Override
  public int getLengthOnReference() {
    return mLengthOnReference;
  }

  @Override
  public int getStartPosition() {
    return mStartPosition;
  }

  @Override
  public String getReferenceName() {
    return mReferenceName;
  }

  @Override
  public int getReferenceId() {
    return mRefId;
  }

  @Override
  public boolean hasReference() {
    return true;
  }

  @Override
  public int getBinNum() {
    return mBin;
  }

  @Override
  public boolean hasCoordinates() {
    return true;
  }

  @Override
  public boolean isUnmapped() {
    return false;
  }

  @Override
  public List<String> getSequenceNames() {
    return mNamesList;
  }

  @Override
  public void close() throws IOException {
    mReader.close();
  }

  @Override
  public boolean hasNext() {
    return mNextLine != null;
  }

  @Override
  public long getVirtualOffset() {
    return mCurrentVirtualOffset;
  }

  @Override
  public long getNextVirtualOffset() {
    return mNextVirtualOffset;
  }

  /**
   * Get the value from the given column of the current record
   * @param col zero base column
   * @return the value
   */
  protected String getColumn(int col) throws IOException {
    populateTabs(col);
    final int pos1 = mTabs[col];
    final int pos2 = mTabs[col + 1];
    if (pos1 >= mCurrentLine.length()) {
      throw new IOException("Data file did not contain expected column " + (col + 1) + " on line: " + mCurrentLine);
    }
    return mCurrentLine.substring(pos1 + 1, pos2);
  }

  private void populateTabs(int colNo) {
    while (mTabsUsed <= colNo + 1) {
      int pos = mTabs[mTabsUsed - 1] + 1;
      while (pos < mCurrentLine.length() && mCurrentLine.charAt(pos) != '\t') {
        pos++;
      }
      mTabs[mTabsUsed++] = pos;
    }
  }

  private void populateNext() {
    do {
      mNextVirtualOffset = mReader.getFilePointer();
      mNextLine = mReader.readLine();
    } while (mNextLine != null && (mNextLine.length() == 0 || mNextLine.charAt(0) == mMeta));
  }
}
