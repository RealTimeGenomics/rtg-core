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


  MatchResult(int size) {
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
    ++mCount;
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
    for (int i = 0; i < size(); ++i) {
      sb.append("templateId: ").append(mTemplateId.getInt(i)).append(" position: ").append(mPosition.getInt(i)).append(" encodedReadId: ").append(mEncodedReadId.getInt(i)).append(" reverse: ").append(mReverse.get(i) > 0).append(StringUtils.LS);
    }
    return sb.toString();
  }
  private static class MatchResultSortProxy implements QuickSort.SortProxy {

    private final MatchResult mResult;
    MatchResultSortProxy(MatchResult result) {
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
