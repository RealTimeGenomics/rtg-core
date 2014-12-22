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
package com.rtg.ngs;

import com.rtg.util.QuickSort;
import com.rtg.util.StringUtils;
import com.rtg.util.array.bitindex.BitCreate;
import com.rtg.util.array.bitindex.BitIndex;
import com.rtg.util.array.intindex.IntChunks;


/**
 * Holds some information about a match result
 */
final class MatchResult {
  IntChunks mTemplateId;
  IntChunks mPosition;
  IntChunks mEncodedReadId;
  BitIndex mReverse;

  long mCount = 0;


  public MatchResult(int size) {
    final int newSize = Math.max(2, size);
    mTemplateId = new IntChunks(newSize);
    mPosition = new IntChunks(newSize);
    mEncodedReadId = new IntChunks(newSize);
    mReverse = BitCreate.createIndex(newSize, 1);
    assert mTemplateId.integrity();
  }

  // Remove references to all the memory consuming data structures
  void evacuateTheBuilding() {
    mTemplateId = null;
    mPosition = null;
    mEncodedReadId = null;
    mReverse = null;
  }

  /**
   * Constructor
   * @param encodedReadId as used by paired end system
   * @param reverse true if reverse complement
   * @param templateId id of template
   * @param position position on template
   */
  public void addMatchResult(int templateId, int position, int encodedReadId, boolean reverse) {
    if (templateId < 0 || encodedReadId < 0) {
      throw new IllegalArgumentException();
    }
    ensureCapacity();
    mEncodedReadId.setInt(mCount, encodedReadId);
    mReverse.set(mCount, reverse ? 1 : 0);
    mTemplateId.setInt(mCount, templateId);
    mPosition.setInt(mCount, position);
    mCount++;
  }

  private void ensureCapacity() {
    if (mCount == mTemplateId.length()) {
      final long extend = mCount / 2;
      mTemplateId.extendBy(extend);
      mPosition.extendBy(extend);
      mEncodedReadId.extendBy(extend);
      mReverse.extendBy(extend);
    }
  }

  public int getTemplateId(long index) {
    return mTemplateId.getInt(index);
  }

  public int getEncodedReadId(long index) {
    return mEncodedReadId.getInt(index);
  }

  public int getPosition(long index) {
    return mPosition.getInt(index);
  }

  public boolean isReverse(long index) {
    return mReverse.get(index) > 0;
  }

  public long size() {
    return mCount;
  }

  public int compare(int i1, int i2) {
    if (mTemplateId.getInt(i1) != mTemplateId.getInt(i2)) {
      return mTemplateId.getInt(i1) < mTemplateId.getInt(i2) ? -1 : 1;
    }
    if (mPosition.getInt(i1) != mPosition.getInt(i2)) {
      return mPosition.getInt(i1) < mPosition.getInt(i2) ? -1 : 1;
    }
    if (mEncodedReadId.getInt(i1) != mEncodedReadId.getInt(i2)) {
      return mEncodedReadId.getInt(i1) < mEncodedReadId.getInt(i2) ? -1 : 1;
    }
    return 0;
  }

  public void swap(int i1, int i2) {
    mTemplateId.swap(i1, i2);
    mPosition.swap(i1, i2);
    mEncodedReadId.swap(i1, i2);
    mReverse.swap(i1, i2);
  }

  public void sort() {
    QuickSort.sort(new MatchResultSortProxy(this));
  }

  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < size(); i++) {
      sb.append("templateId: ").append(mTemplateId.getInt(i)).append(" position: ").append(mPosition.getInt(i)).append(" encodedReadId: ").append(mEncodedReadId.getInt(i)).append(" reverse: ").append(mReverse.get(i) > 0).append(StringUtils.LS);
    }
    return sb.toString();
  }
  private static class MatchResultSortProxy implements QuickSort.SortProxy {

    private final MatchResult mResult;
    public MatchResultSortProxy(MatchResult result) {
      mResult = result;
    }

    @Override
    public int compare(long index1, long index2) {
      return mResult.compare((int) index1, (int) index2);
    }

    @Override
    public long length() {
      return mResult.size();
    }

    @Override
    public void swap(long index1, long index2) {
      mResult.swap((int) index1, (int) index2);
    }

  }

}
