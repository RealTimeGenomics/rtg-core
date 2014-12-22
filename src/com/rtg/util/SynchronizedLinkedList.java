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

import java.util.Collection;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * A synchronized linked list wrapper, implementing <code>Deque</code>
 * @param <E> the type of elements held in this collection
 */
public class SynchronizedLinkedList<E> implements Deque<E> {

  private final LinkedList<E> mList;
  private final Object mMutex;

  /**
   * Create a synchronized linked list
   * @param list the linked list to wrap
   */
  public SynchronizedLinkedList(LinkedList<E> list) {
    mList = list;
    mMutex = this;
  }

  @Override
  public boolean isEmpty() {
    synchronized (mMutex) {
      return mList.isEmpty();
    }
  }

  @Override
  public Object[] toArray() {
    synchronized (mMutex) {
      return mList.toArray();
    }
  }

  @Override
  public <T> T[] toArray(T[] a) {
    synchronized (mMutex) {
      return mList.toArray(a);
    }
  }

  @Override
  public boolean containsAll(Collection<?> c) {
    synchronized (mMutex) {
      return mList.containsAll(c);
    }
  }

  @Override
  public boolean addAll(Collection<? extends E> c) {
    synchronized (mMutex) {
      return mList.addAll(c);
    }
  }

  @Override
  public boolean removeAll(Collection<?> c) {
    synchronized (mMutex) {
      return mList.removeAll(c);
    }
  }

  @Override
  public boolean retainAll(Collection<?> c) {
    synchronized (mMutex) {
      return mList.retainAll(c);
    }
  }

  @Override
  public void clear() {
    synchronized (mMutex) {
      mList.clear();
    }
  }

  @Override
  public void addFirst(E e) {
    synchronized (mMutex) {
      mList.addFirst(e);
    }
  }

  @Override
  public void addLast(E e) {
    synchronized (mMutex) {
      mList.addLast(e);
    }
  }

  @Override
  public boolean offerFirst(E e) {
    synchronized (mMutex) {
      return mList.offerFirst(e);
    }
  }

  @Override
  public boolean offerLast(E e) {
    synchronized (mMutex) {
      return mList.offerLast(e);
    }
  }

  @Override
  public E removeFirst() {
    synchronized (mMutex) {
      return mList.removeFirst();
    }
  }

  @Override
  public E removeLast() {
    synchronized (mMutex) {
      return mList.removeLast();
    }
  }

  @Override
  public E pollFirst() {
    synchronized (mMutex) {
      return mList.pollFirst();
    }
  }

  @Override
  public E pollLast() {
    synchronized (mMutex) {
      return mList.pollLast();
    }
  }

  @Override
  public E getFirst() {
    synchronized (mMutex) {
      return mList.getFirst();
    }
  }

  @Override
  public E getLast() {
    synchronized (mMutex) {
      return mList.getLast();
    }
  }

  @Override
  public E peekFirst() {
    synchronized (mMutex) {
      return mList.peekFirst();
    }
  }

  @Override
  public E peekLast() {
    synchronized (mMutex) {
      return mList.peekLast();
    }
  }

  @Override
  public boolean removeFirstOccurrence(Object o) {
    synchronized (mMutex) {
      return mList.removeFirstOccurrence(o);
    }
  }

  @Override
  public boolean removeLastOccurrence(Object o) {
    synchronized (mMutex) {
      return mList.removeLastOccurrence(o);
    }
  }

  @Override
  public boolean add(E e) {
    synchronized (mMutex) {
      return mList.add(e);
    }
  }

  @Override
  public boolean offer(E e) {
    synchronized (mMutex) {
      return mList.offer(e);
    }
  }

  @Override
  public E remove() {
    synchronized (mMutex) {
      return mList.remove();
    }
  }

  @Override
  public E poll() {
    synchronized (mMutex) {
      return mList.poll();
    }
  }

  @Override
  public E element() {
    synchronized (mMutex) {
      return mList.element();
    }
  }

  @Override
  public E peek() {
    synchronized (mMutex) {
      return mList.peek();
    }
  }

  @Override
  public void push(E e) {
    synchronized (mMutex) {
      mList.push(e);
    }
  }

  @Override
  public E pop() {
    synchronized (mMutex) {
      return mList.pop();
    }
  }

  @Override
  public boolean remove(Object o) {
    synchronized (mMutex) {
      return mList.remove(o);
    }
  }

  @Override
  public boolean contains(Object o) {
    synchronized (mMutex) {
      return mList.contains(o);
    }
  }

  @Override
  public int size() {
    synchronized (mMutex) {
      return mList.size();
    }
  }

  @Override
  public Iterator<E> iterator() {
    synchronized (mMutex) {
      return mList.iterator();
    }
  }

  @Override
  public Iterator<E> descendingIterator() {
    synchronized (mMutex) {
      return mList.descendingIterator();
    }
  }

  @Override
  public String toString() {
    synchronized (mMutex) {
      return mList.toString();
    }
  }
}
