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

package com.rtg.assembler;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.iterators.IteratorHelper;

/**
 */
public class HashMapDeBruijnGraph implements Iterable<Kmer>, DeBruijnGraph {
  /**
   * This map treats all Kmer values as being equivalent to their reverse complement
   * Override a couple of the methods to enforce this.
   * Am I being too clever
   * Should probably come back and encapsulate this to hide the methods that haven't been overridden correctly.
   */
  final HashMap<Kmer, DeBruijnNode> mDeBruijnGraph =  new ComplementHashMap();
  int mThreshold = 0;

  @Override
  public long bytes() {
    return 0; // Size estimate not available
  }

  private static class ComplementHashMap extends HashMap<Kmer, DeBruijnNode> {
    @Override
    public DeBruijnNode put(Kmer key, DeBruijnNode value) {
      return super.put(key.minimalKmer(), value);
    }
    @Override
    public DeBruijnNode get(Object key) {
      if (!(key instanceof Kmer)) {
        throw new IllegalArgumentException();
      }
      final Kmer kmer = (Kmer) key;
      final Kmer minimal = kmer.minimalKmer();
      return super.get(minimal);
    }
    @Override
    public boolean containsKey(Object key) {
      if (!(key instanceof Kmer)) {
        throw new IllegalArgumentException();
      }
      final Kmer kmer = (Kmer) key;
      return super.containsKey(kmer.minimalKmer());
    }

  }

  HashMapDeBruijnGraph(KmerIterableFactoryInterface factory) {
    try (final KmerIterable kmers = factory.makeIterable()) {
      for (final Kmer k : kmers) {
        add(k);
      }
    } catch (IOException e) {
      throw new NoTalkbackSlimException(e, ErrorType.IO_ERROR);
    }
  }


  private void add(Kmer kmer) {
    if (mDeBruijnGraph.containsKey(kmer)) {
      final DeBruijnNode node = mDeBruijnGraph.get(kmer);
      node.mCount++;
    } else {
      mDeBruijnGraph.put(kmer, new DeBruijnNode());
    }
  }

  @Override
  public int frequency(Kmer k) {
    return getNode(k).mCount;
  }

  DeBruijnNode getNode(Kmer key) {
    return mDeBruijnGraph.get(key);
  }

  /**
   * set the threshold used to determine if a hash is good
   * @param goodThreshold number of times a hash should occur before you'll believe it comes from the genome rather than
   *                      read error
   */
  @Override
  public void setThreshold(int goodThreshold) {
    mThreshold = goodThreshold;
  }
  @Override
  public void setBuilt(Kmer k, boolean built) {
    getNode(k).mPartOfContig = built;
  }
  @Override
  public boolean isBuilt(Kmer k) {
    return getNode(k).mPartOfContig;
  }
  @Override
  public boolean contains(Kmer k) {
    return mDeBruijnGraph.containsKey(k) && frequency(k) > mThreshold;
  }

  /**
   * The values stored in our Kmer map
   * Each node in the graph has a count and may one day belong to a <code>PreContig</code>
   */
  public static class DeBruijnNode {
    int mCount = 1;
    boolean mPartOfContig = false;

    @Override
    public String toString() {
      return "partOfContig=" + mPartOfContig + " mCount=" + mCount;
    }
  }

  @Override
  public Iterator<Kmer> iterator() {
    final Set<Entry<Kmer, DeBruijnNode>> set = mDeBruijnGraph.entrySet();
    final Iterator<Entry<Kmer, DeBruijnNode>> it = set.iterator();
    return new KmerIterator(it, mThreshold);
  }

  private static class KmerIterator extends IteratorHelper<Kmer> {
    private Map.Entry<Kmer, DeBruijnNode> mNext = null;
    private final Iterator<Map.Entry<Kmer, DeBruijnNode>> mInternal;
    private final int mThreshold;

    KmerIterator(Iterator<Map.Entry<Kmer, DeBruijnNode>> internalIterator, int threshold) {
      mInternal = internalIterator;
      mThreshold = threshold;
    }

    @Override
    protected void step() {
      mNext = mInternal.hasNext() ? mInternal.next() : null;
    }

    @Override
    protected boolean atEnd() {
      return mNext == null && !mInternal.hasNext();
    }

    @Override
    protected boolean isOK() {
      return mNext != null && mNext.getValue().mCount > mThreshold;
    }

    @Override
    protected Kmer current() {
      return mNext.getKey();
    }
  }

}
