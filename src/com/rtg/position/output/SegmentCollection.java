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
    mSize++;
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
    for (int i = 0; i < mSize; i++) {
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
    for (int i = 0; i < mSize; i++) {
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
    mSize--;
    final Segment res = mSegments[mSize];
    mSegments[mSize] = null;
    return res;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SegmentCollection [").append(mSize).append("]").append(StringUtils.LS);
    for (int i = 0; i < mSize; i++) {
      sb.append("[").append(i).append("] ").append(mSegments[i].toString()).append(StringUtils.LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mSize; i++) {
      Exam.assertTrue(mSegments[i] != null);
      if (i > 0) {
        Exam.assertEquals(mSegments[i].isEmpty(), mSegments[i - 1].isEmpty());
      }
    }
    for (int i = mSize; i < mSegments.length; i++) {
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
