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
package com.rtg.reader;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import com.rtg.mode.SequenceType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.InformationType;

/**
 * Concatenates multiple <code>SequenceDataSource</code>s
 *
 * @param <T> <code>SequenceDataSource</code> or a subclass
 */
public class ConcatSequenceDataSource<T extends SequenceDataSource> implements SequenceDataSource {

  private final List<T> mSources;
  private final Iterator<T> mIterator;
  private final List<String> mNames;
  private SequenceDataSource mCurrent;
  private long mWarningCount;
  private long mDustedCount;
  private int mNameIndex;
  private long mMinLength = Long.MAX_VALUE;
  private long mMaxLength = Long.MIN_VALUE;

  /**
   *
   * @param sources <code>SequencesDataSource</code>s to read
   * @param names names of sources to print to user (pass <code>null</code> to not print)
   */
  public ConcatSequenceDataSource(final List<T> sources, final List<String> names) {
    if (sources == null || sources.size() == 0) {
      throw new IllegalArgumentException("Cannot concatenate 0 sources");
    }
    mSources = sources;
    mIterator = mSources.iterator();
    mNames = names;
    mCurrent = mIterator.next();
    mWarningCount = 0;
    mNameIndex = 0;
    mDustedCount = 0;
    printFileNumber();
  }

  @Override
  public SequenceType type() {
    return mCurrent.type();
  }

  @Override
  public boolean nextSequence() throws IOException {
    while (!mCurrent.nextSequence()) {
      if (!mIterator.hasNext()) {
        return false;
      }
      mWarningCount += mCurrent.getWarningCount();
      mDustedCount += mCurrent.getDusted();
      mMinLength = Math.min(mMinLength, mCurrent.getMinLength());
      mMaxLength = Math.max(mMaxLength, mCurrent.getMaxLength());
      mCurrent.close();
      mCurrent = mIterator.next();
      mNameIndex++;
      printFileNumber();
    }
    return true;
  }

  @Override
  public String name() throws IOException {
    return mCurrent.name();
  }

  @Override
  public byte[] sequenceData() throws IOException {
    return mCurrent.sequenceData();
  }

  @Override
  public byte[] qualityData() throws IOException {
    return mCurrent.qualityData();
  }

  @Override
  public boolean hasQualityData() {
    return mCurrent.hasQualityData();
  }

  @Override
  public int currentLength() throws IOException {
    return mCurrent.currentLength();
  }

  @Override
  public void close() throws IOException {
    mCurrent.close();
    while (mIterator.hasNext()) {
      mIterator.next().close();
    }
  }

  @Override
  public void setDusting(boolean val) {
    for (final SequenceDataSource source : mSources) {
      source.setDusting(val);
    }
  }

  @Override
  public long getWarningCount() {
    return mWarningCount + mCurrent.getWarningCount();
  }

  private void printFileNumber() {
    if (mNames != null) {
      Diagnostic.info(InformationType.PROCESSING_ITEM_N_OF_N, true, "", mNames.get(mNameIndex), Integer.toString(mNameIndex + 1), Integer.toString(mSources.size()));
    }
  }

  @Override
  public long getDusted() {
    return mDustedCount + mCurrent.getDusted();
  }

  @Override
  public long getMaxLength() {
    return Math.max(mMaxLength, mCurrent.getMaxLength());
  }

  @Override
  public long getMinLength() {
    return Math.min(mMinLength, mCurrent.getMinLength());
  }
}
