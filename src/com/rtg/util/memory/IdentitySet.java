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
package com.rtg.util.memory;

import java.util.NoSuchElementException;

import com.rtg.util.diagnostic.Diagnostic;


/**
 * A set implementation (it implements a subset of java.util.Set). It
 * remembers all objects using identity (==) for comparison and hashes
 * on the default <code>hashCode()</code> method. This is in distinction to the
 * java.util.Set implementations which use equals for comparison and
 * use the objects own hash function.
 *
 */
public class IdentitySet {

  /** Default number of initial entries. */
  static final int INITIAL_SIZE = 101; // a prime roughly the right size

  /** Minimum length permitted (stops stupid stuff at small values). */
  static final int MIN_SIZE = 11;

  /**
   * If more than this fraction of entries are added then the system is
   * resized.
   */
  static final double LOAD_FACTOR = 0.3;

  /**
   * Number of entries that are searched (in sequence) before
   * overflowing.
   */
  static final int BLOCK_SIZE = 8; //handy power of 2

  /** Array with pointers to objects. */
  private Object[] mObj;

  /** Nominal length of <code>mObj</code>. */
  private int mLength;

  /**
   * True length of <code>mObj</code> (rounded up to allow for searching without
   * having to check if off end of array)
   */
  private int mTrueLength;

  /**
   * A recursive IdentitySet used to handle overflows (that is entries
   * not encountered before searching BLOCK_SIZE entries.
   */
  private IdentitySet mOverflow = null;

  /** Number of objects stored. */
  private int mN = 0;

  static final int FOUND = 0;
  static final int EMPTY = 1;
  static final int OVER_FLOW = 2;

  /** If R_STATE == FOUND or EMPTY then the relevant index position. */
  private int mRindex;

  /**
   * @param length initial number of entries in <code>mObj</code>.
   */
  public IdentitySet(final int length) {
    if (length < MIN_SIZE) {
      mLength = MIN_SIZE;
    } else {
      mLength = length;
    }
    mTrueLength = mLength + BLOCK_SIZE - 1;
    mObj = new Object[mTrueLength];
    assert integrity();
  }


  IdentitySet() {
    this(INITIAL_SIZE);
  }


  /**
   * Adds all items from an array to the set.
   *
   * @param objs an <code>Object[]</code> value
   */
  public void addAll(final Object[] objs) {
    for (Object obj : objs) {
      if (obj != null) {
        add(obj);
      }
    }
  }

  /**
   * Add a new object to the set.
   *
   * @param obj Object to be added.
   * @return true iff <code>obj</code> not already present.
   */
  public boolean add(final Object obj) {
    if (obj == null) {
      throw new IllegalArgumentException();
    }

    final int state = find(obj);
    final boolean ans;
    switch (state) {
        case FOUND:
          ans = false;
          break;
        case EMPTY:
          try {
            if (mObj[mRindex] != null) {
              throw new RuntimeException("Borken");
            }
            mObj[mRindex] = obj;
          } catch (final Exception e) {
            Diagnostic.userLog(">>> " + obj.getClass().getName());
            Diagnostic.userLog(e);
          }
          ++mN;
          if (mN > mLength * LOAD_FACTOR) {
            resize();
          }
          ans = true;
          break;
        case OVER_FLOW: //should be very rare
          if (mOverflow == null) {
            mOverflow = new IdentitySet();
          }
          ans = mOverflow.add(obj);
          break;
        default:
          throw new RuntimeException();
    }
    assert integrity();
    return ans;
  }


  private void resize() {
    final IdentitySet newset = new IdentitySet(mLength * 2 + 1);
    for (final java.util.Iterator<Object> iter = getIterator(); iter.hasNext(); ) {
      final Object obj = iter.next();
      newset.add(obj);
    }
    assert newset.integrity();
    mObj = newset.mObj;
    mLength = newset.mLength;
    mTrueLength = newset.mTrueLength;
    mOverflow = newset.mOverflow;
    mN = newset.mN;
    assert integrity();
  }


  /**
   * Search for an object. Results placed in r_ variables.
   *
   * @param obj the object being sought.
   * @return state
   */
  private int find(final Object obj) {
    final int prb = probe(obj);
    return find(obj, prb);
  }


  /*
   * Search for an object. Results placed in r_ variables.
   */
  private int find(final Object obj, final int prb) {
    mRindex = prb % mLength;
    for (int i = 0; i < BLOCK_SIZE; ++i, ++mRindex) {
      final Object vobj = mObj[mRindex];
      if (vobj == null) {
        //sRstate = EMPTY;
        return EMPTY;
      }
      if (vobj == obj) {
        //sRstate = FOUND;
        return FOUND;
      }
    }

    if (mOverflow == null) {
      //sRstate = OVER_FLOW;
    } else {
      if (mOverflow.find(obj, prb) == FOUND) {
        return FOUND;
      }
    }
    return OVER_FLOW;
  }


  private int probe(final Object obj) {
    final int x = System.identityHashCode(obj) * 11111111;
    return x < 0 ? -x : x;
  }


  /**
   * Membership test.
   *
   * @param obj object to test
   * @return true if a member
   */
  public boolean contains(final Object obj) {
    return find(obj) == FOUND;
  }


  /**
   * Number of entries in set.
   *
   * @return Number of entries in set.
   */
  public int size() {
    return mN + (mOverflow == null ? 0 : mOverflow.size());
  }


  /**
   * Check if any entries are present.
   *
   * @return true iff there are no entries present.
   */
  public boolean isEmpty() {
    return mN == 0;
  }


  /**
   * Used for implementing iterator() method.
   *
   */
  private static class Iterator implements java.util.Iterator<Object> {
    private IdentitySet mSet;
    private int mIndex = 0;


    Iterator(final IdentitySet set) {
      mSet = set;
      check();
    }

    @Override
    public boolean hasNext() {
      return mSet != null;
    }


    @Override
    public Object next() {
      final Object result = mSet.mObj[mIndex];
      if (result == null) {
        throw new NoSuchElementException();
      }
      ++mIndex;
      check();
      return result;
    }


    private void check() {
      while (true) {
        if (mIndex == mSet.mTrueLength) {
          mIndex = 0;
          mSet = mSet.mOverflow;
          if (mSet == null) {
            return;
          } else {
            continue;
          }
        } else if (mSet.mObj[mIndex] != null) {
          return;
        } else {
          ++mIndex;
          continue;
        }
      }
    } //check


    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }


  /**
   * Iterator over all entries in set.
   *
   * @return Iterator over all entries in set.
   */
  public java.util.Iterator<Object> getIterator() {
    return new Iterator(this);
  }


  boolean integrity() {
    assert INITIAL_SIZE >= MIN_SIZE;
    assert LOAD_FACTOR > 0.0 && LOAD_FACTOR < 1.0;
    assert BLOCK_SIZE > 1;
    assert mObj != null;
    assert mObj.length == mTrueLength;
    assert mTrueLength == (mLength + BLOCK_SIZE - 1);
    assert mN >= 0;

    //chekc that nothing more than BLOCK_SIZE away from initial probe point
    int lastNull = -1;
    for (int i = 0; i < mObj.length; ++i) {
      if (mObj[i] == null) {
        lastNull = i;
      } else {
        final int p = probe(mObj[i]) % mLength;
        assert p > lastNull && p <= i && i - BLOCK_SIZE <= p : i + ":" + p + ":" + lastNull;
      }
    }
    if (mOverflow != null) {
      mOverflow.integrity();
      //check if all overflowed objects really were supposed to overflow
      for (final java.util.Iterator<Object> iter = mOverflow.getIterator(); iter.hasNext(); ) {
        final Object obj = iter.next();
        final int p = probe(obj) % mLength;
        for (int i = p; i < p + BLOCK_SIZE; ++i) {
          assert mObj[i] != null;
        }
      }
    }
    int count = 0;
    for (final java.util.Iterator<Object> iter = getIterator(); iter.hasNext(); ) {
      final Object obj = iter.next();
      assert obj != null;
      ++count;
    }
    assert size() == count;
    return true;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    toString(sb);
    return sb.toString();
  }


  void toString(final StringBuilder sb) {
    sb.append("IdentitySet:\n");
    for (final java.util.Iterator<Object> iter = getIterator(); iter.hasNext(); ) {
      final Object obj = iter.next();
      sb.append(obj.toString());
      sb.append("\n");
    }
    sb.append("End IdentitySet\n");
  }



}

