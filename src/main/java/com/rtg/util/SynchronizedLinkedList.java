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

  private final LinkedList<E> mLinkedList;
  private final Object mMutex;

  /**
   * Create a synchronized linked list
   * @param linkedList the linked list to wrap
   */
  public SynchronizedLinkedList(LinkedList<E> linkedList) {
    mLinkedList = linkedList;
    mMutex = this;
  }

  @Override
  public boolean isEmpty() {
    synchronized (mMutex) {
      return mLinkedList.isEmpty();
    }
  }

  @Override
  public Object[] toArray() {
    synchronized (mMutex) {
      return mLinkedList.toArray();
    }
  }

  @Override
  public <T> T[] toArray(T[] a) {
    synchronized (mMutex) {
      return mLinkedList.toArray(a);
    }
  }

  @Override
  public boolean containsAll(Collection<?> c) {
    synchronized (mMutex) {
      return mLinkedList.containsAll(c);
    }
  }

  @Override
  public boolean addAll(Collection<? extends E> c) {
    synchronized (mMutex) {
      return mLinkedList.addAll(c);
    }
  }

  @Override
  public boolean removeAll(Collection<?> c) {
    synchronized (mMutex) {
      return mLinkedList.removeAll(c);
    }
  }

  @Override
  public boolean retainAll(Collection<?> c) {
    synchronized (mMutex) {
      return mLinkedList.retainAll(c);
    }
  }

  @Override
  public void clear() {
    synchronized (mMutex) {
      mLinkedList.clear();
    }
  }

  @Override
  public void addFirst(E e) {
    synchronized (mMutex) {
      mLinkedList.addFirst(e);
    }
  }

  @Override
  public void addLast(E e) {
    synchronized (mMutex) {
      mLinkedList.addLast(e);
    }
  }

  @Override
  public boolean offerFirst(E e) {
    synchronized (mMutex) {
      return mLinkedList.offerFirst(e);
    }
  }

  @Override
  public boolean offerLast(E e) {
    synchronized (mMutex) {
      return mLinkedList.offerLast(e);
    }
  }

  @Override
  public E removeFirst() {
    synchronized (mMutex) {
      return mLinkedList.removeFirst();
    }
  }

  @Override
  public E removeLast() {
    synchronized (mMutex) {
      return mLinkedList.removeLast();
    }
  }

  @Override
  public E pollFirst() {
    synchronized (mMutex) {
      return mLinkedList.pollFirst();
    }
  }

  @Override
  public E pollLast() {
    synchronized (mMutex) {
      return mLinkedList.pollLast();
    }
  }

  @Override
  public E getFirst() {
    synchronized (mMutex) {
      return mLinkedList.getFirst();
    }
  }

  @Override
  public E getLast() {
    synchronized (mMutex) {
      return mLinkedList.getLast();
    }
  }

  @Override
  public E peekFirst() {
    synchronized (mMutex) {
      return mLinkedList.peekFirst();
    }
  }

  @Override
  public E peekLast() {
    synchronized (mMutex) {
      return mLinkedList.peekLast();
    }
  }

  @Override
  public boolean removeFirstOccurrence(Object o) {
    synchronized (mMutex) {
      return mLinkedList.removeFirstOccurrence(o);
    }
  }

  @Override
  public boolean removeLastOccurrence(Object o) {
    synchronized (mMutex) {
      return mLinkedList.removeLastOccurrence(o);
    }
  }

  @Override
  public boolean add(E e) {
    synchronized (mMutex) {
      return mLinkedList.add(e);
    }
  }

  @Override
  public boolean offer(E e) {
    synchronized (mMutex) {
      return mLinkedList.offer(e);
    }
  }

  @Override
  public E remove() {
    synchronized (mMutex) {
      return mLinkedList.remove();
    }
  }

  @Override
  public E poll() {
    synchronized (mMutex) {
      return mLinkedList.poll();
    }
  }

  @Override
  public E element() {
    synchronized (mMutex) {
      return mLinkedList.element();
    }
  }

  @Override
  public E peek() {
    synchronized (mMutex) {
      return mLinkedList.peek();
    }
  }

  @Override
  public void push(E e) {
    synchronized (mMutex) {
      mLinkedList.push(e);
    }
  }

  @Override
  public E pop() {
    synchronized (mMutex) {
      return mLinkedList.pop();
    }
  }

  @Override
  public boolean remove(Object o) {
    synchronized (mMutex) {
      return mLinkedList.remove(o);
    }
  }

  @Override
  public boolean contains(Object o) {
    synchronized (mMutex) {
      return mLinkedList.contains(o);
    }
  }

  @Override
  public int size() {
    synchronized (mMutex) {
      return mLinkedList.size();
    }
  }

  @Override
  public Iterator<E> iterator() {
    synchronized (mMutex) {
      return mLinkedList.iterator();
    }
  }

  @Override
  public Iterator<E> descendingIterator() {
    synchronized (mMutex) {
      return mLinkedList.descendingIterator();
    }
  }

  @Override
  public String toString() {
    synchronized (mMutex) {
      return mLinkedList.toString();
    }
  }
}
