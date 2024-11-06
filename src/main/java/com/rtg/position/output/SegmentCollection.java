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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;



/**
 * Holds a collection of <code>Segment</code>s.
 * They can be added and removed one by one.
 */
public class SegmentCollection extends IntegralAbstract {

  private final Segment[] mSegments;

  private int mSize = 0;

  /**
   * @param size maximum number of <code>Segment</code>s in the collection.
   */
  public SegmentCollection(final int size) {
    mSegments = new Segment[size];
  }

  /**
   * Get the number of entries in the collection.
   * @return the number of entries in the collection.
   */
  public int size() {
    return mSize;
  }

  /**
   * Add a segment.
   * @param segment It is not allowed to be null.
   */
  public void add(final Segment segment) {
    assert segment != null;
//    assert mSize == 0 || segment.compareTo(mSegments[mSize - 1]) > 0 || segment.isEmpty() : segment + ":" + mSegments[mSize - 1];
    mSegments[mSize] = segment;
    ++mSize;
  }

  /**
   * Get the <code>Segment</code> specified by index.
   * @param index specifies the <code>Segment</code>.
   * @return the specified <code>Segment</code> (can never be null).
   */
  public Segment get(final int index) {
    if (index >= mSize || index < 0) {
      throw new IndexOutOfBoundsException(index + ":" + mSize);
    }
    return mSegments[index];
  }

  /**
   * Remove all <code>Segment</code>s.
   */
  public void clear() {
    for (int i = 0; i < mSize; ++i) {
      mSegments[i] = null;
    }
    mSize = 0;
  }

  /**
   * Write all the segments also freeing them and clearing collection.
   * @param out where to write the output.
   * @param free list of free <code>Segment</code>s.
   * @param searchPosition position in the search sequence.
   * @throws IOException if error writing output.
   */
  public void flush(final SegmentWriter out, final SegmentCollection free, final int searchPosition) throws IOException {
    for (int i = 0; i < mSize; ++i) {
      final Segment seg = mSegments[i];
      out.write(seg, searchPosition);
      seg.clear();
      free.add(seg);
      mSegments[i] = null;
    }
    mSize = 0;
  }

  /**
   * Remove a segment from end of collection and return it.
   * Useful for free list.
   * @return removed segment.
   */
  public Segment removeNext() {
    if (mSize == 0) {
      return null;
    }
    --mSize;
    final Segment res = mSegments[mSize];
    mSegments[mSize] = null;
    return res;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SegmentCollection [").append(mSize).append("]").append(StringUtils.LS);
    for (int i = 0; i < mSize; ++i) {
      sb.append("[").append(i).append("] ").append(mSegments[i]).append(StringUtils.LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mSize; ++i) {
      Exam.assertTrue(mSegments[i] != null);
      if (i > 0) {
        Exam.assertEquals(mSegments[i].isEmpty(), mSegments[i - 1].isEmpty());
      }
    }
    for (int i = mSize; i < mSegments.length; ++i) {
      Exam.assertTrue(mSegments[i] == null);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSize >= 0 && mSize <= mSegments.length);
    return true;
  }

}
