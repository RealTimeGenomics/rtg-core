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
import com.rtg.util.memory.MemoryUsage;

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
    return MemoryUsage.size(mDeBruijnGraph);
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
