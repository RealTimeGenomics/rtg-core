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

package com.rtg.util;

import java.util.Iterator;

/**
 * @param <T> the type of object to collect
 */
public class BasicLinkedListNode<T> implements Iterable<T> {
  final T mValue;
  final BasicLinkedListNode<T> mNext;
  final int mSize;

  /**
   * Create a new pointer to the head of a list
   * @param value the value to store at the head
   * @param tail the rest of the list
   */
  public BasicLinkedListNode(T value, BasicLinkedListNode<T> tail) {
    mNext = tail;
    mValue = value;
    if (tail == null) {
      mSize = 1;
    } else {
      mSize = tail.mSize + 1;
    }
  }

  @Override
  public Iterator<T> iterator() {
    return new BasicLinkedListNodeIterator<>(this);
  }

  /**
   * @return the number of elements in this list
   */
  public int size() {
    return mSize;
  }

  /**
   * @return the value stored in the node
   */
  public T getValue() {
    return mValue;
  }
}
