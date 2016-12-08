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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Path;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Histogram;
import com.rtg.util.Pair;
import com.rtg.util.array.CommonIndex;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;

/**
 * Construct a <code>DeBruijn</code> graph of Kmer nodes
 */
public class DeBruijnGraphBuilder {

  private static final boolean FORCE_HASHMAP = false; //Boolean.valueOf(System.getProperty("rtg.assembler.usehashmap", "false"));

  /**
   * This map treats all Kmer values as being equivalent to their reverse complement
   * Override a couple of the methods to enforce this.
   * Am I being too clever
   * Should probably come back and encapsulate this to hide the methods that haven't been overridden correctly.
   */
  final DeBruijnGraph mDeBruijnGraph;

  /**
   * Kmer are required to have at least this many occurrences before we'll believe in them
   * @param goodThreshold kmer with frequency below this value will be ignored
   */
  public void setGoodThreshold(int goodThreshold) {
    mDeBruijnGraph.setThreshold(goodThreshold);
  }


  final int mKmerSize;
  final KmerFactory mFactory;

  /** Padding on length of tips */
  final int mTipConstant;
  /** Map from <code>PreContig.mId</code> to <code>PreContig</code> When we build a list of <code>PreContigs</code> we'll store it here */
  private final GraphKmerAttribute mContigGraph;


  /**
   * Construct a <code>DeBruijn</code> graph with the default Kmer factory
   * @param files iterate over these to get all the reads, <code>byte[]</code>s and then k-mers.
   * @param kmerSize length of Kmer
   * @param tipConstant length beyond 2 kmer that will still be considered a tip
   * @param numberThreads number of threads to run with
   */
  public DeBruijnGraphBuilder(List<ReadPairSource> files, int kmerSize, int tipConstant, int numberThreads) {
    this(files, kmerSize, StringKmer.factory(), tipConstant, numberThreads);
  }
  /**
   * Construct a graph for the given Kmer size
   * @param sources iterate over these to get all the reads, <code>byte[]</code>s and then k-mers.
   * @param kmerSize size of Kmer in bases
   * @param factory a factory that will construct the Kmer objects
   * @param tipConstant length beyond 2 kmer that will still be considered a tip
   * @param numberThreads number of threads to run with
   */
  public DeBruijnGraphBuilder(List<ReadPairSource> sources, int kmerSize, KmerFactory factory, int tipConstant, int numberThreads) {
    mKmerSize = kmerSize;
    mContigGraph = new GraphKmerAttribute(kmerSize - 1, new HashMap<String, String>(), new HashMap<String, String>());
    mFactory = factory;
    mTipConstant = tipConstant;
    final OneShotTimer init = new OneShotTimer("DeBruijn_build");

    if (kmerSize <= 32 && !FORCE_HASHMAP) {
      final long size = size(sources, kmerSize);
      mDeBruijnGraph =  new LowKDeBruijnGraph(new KmerIterableFactory(sources, factory, kmerSize), size, kmerSize);
    } else {
      mDeBruijnGraph =  new HashMapDeBruijnGraph(new KmerIterableFactory(sources, factory, kmerSize));
    }
    init.stopLog();
  }

  /**
   * @param files containing sequences.
   * @param kmerSize length of a kmer (k).
   * @return total number of kmers in all the sequences.
   */
  static long size(final List<ReadPairSource> files, final int kmerSize) {
    long size = 0;
    for (final ReadPairSource source : files) {
      for (SequencesReader sr : source.mReaders) {
        try {
          final int[] sl = sr.sequenceLengths(0, sr.numberSequences());
          for (final int length : sl) {
            final int excess = length - kmerSize + 1;
            if (excess > 0) {
              size += excess;
            }
          }
        } catch (final IOException e) {
          throw new RuntimeException(e);
        }
      }
    }
    return size;
  }

  static int computeThreshold(final Histogram histogram) {
    Diagnostic.info("Hash frequency histogram: " + histogram.toString());
    for (int i = 2; i < histogram.getLength(); ++i) {
      if (histogram.getValue(i) > histogram.getValue(i - 1)) {
        // Empirically determined on E. faecalis
        return Math.max(1, (i - 1) / 2);
      }
    }
    return 1;
  }

  /**
   * Find the first local minimum in the Kmer frequency histogram to be used as a threshold for a "good" Kmer
   * @return the Kmer frequency threshold
   */
  int calculateGoodThreshold() {
    final Histogram histogram = new Histogram();
    for (final Kmer kmer : mDeBruijnGraph) {
      histogram.increment(Math.min(mDeBruijnGraph.frequency(kmer), 1000));
    }
    return computeThreshold(histogram);
  }

  static String complement(String hash) {
    return DnaUtils.reverseComplement(hash);
  }

  GraphKmerAttribute preContigGraph() {
    return mContigGraph;
  }
  /**
   * Mark <code>DeBruijn</code> nodes and accumulate a set of <code>PreContig</code>
   */
  void buildPreContigs() {
    //final Map<Long, PreContig> preContigMap = new HashMap<>();
    long id = 0;
    for (final Kmer kmer: mDeBruijnGraph) {
      final int frequency = mDeBruijnGraph.frequency(kmer);
      if (!mDeBruijnGraph.isBuilt(kmer)) {
        final PreContig pc = new PreContig(id, kmer, frequency);
        final Set<Kmer> visited = new HashSet<>();
        visited.add(kmer);
        walk(pc, false, kmer, visited);
        walk(pc, true, kmer, visited);
        mDeBruijnGraph.setBuilt(kmer, true);
        final long graphId = mContigGraph.addContig(pc);
        mContigGraph.setKmerFreq(graphId, pc.mKmerCount);
        //preContigMap.put(id, pc);
        ++id;
      }
    }
    final Map<Kmer, Long> contigEnds = new HashMap<>();
    for (long i = 1; i <= mContigGraph.numberContigs(); ++i) {
      final Kmer startKmer = startKmer(mContigGraph.contig(i));
      contigEnds.put(startKmer, i);
      final Kmer endKmer = endKmer(mContigGraph.contig(i));
      contigEnds.put(endKmer, i);
    }
    for (long i = 1; i <= mContigGraph.numberContigs(); ++i) {
      final  Contig contig = mContigGraph.contig(i);
      // We ensure we only add each link once by only doing it for the smaller id out of source & destination
      // The condition is >= for the end link and > for the start link on purpose. This is to prevent cyclical contigs
      // from linking to themselves repeatedly.
      final Kmer startKmer = startKmer(contig);
      final Kmer endKmer = endKmer(contig);
      for (byte b = (byte) DNA.A.ordinal(); b <= DNA.T.ordinal(); ++b) {
        final List<Long> links = getLinks(true, endKmer.successor(b), contigEnds);
        for (final long link : links) {
          if (Math.abs(link) >= i) {
            final Path p = new PathArray(i, (link > 0 ? 1 : -1) * Math.abs(link));
            mContigGraph.addPath(p);
          }
        }
      }

      for (byte b = (byte) DNA.A.ordinal(); b <= DNA.T.ordinal(); ++b) {
        final List<Long> links = getLinks(false, startKmer.predecessor(b), contigEnds);
        for (final long link : links) {
          if (Math.abs(link) > i) {
            final Path p = new PathArray(-i, (link > 0 ? -1 : 1) * Math.abs(link));
            mContigGraph.addPath(p);
          }
        }
      }
    }

    Diagnostic.developerLog("mContigGraph mem size:" + mContigGraph.bytes());
  }

  private Kmer endKmer(Contig contig) {
    return mFactory.make(contig, contig.length() - mKmerSize, contig.length());
  }

  private Kmer startKmer(Contig contig) {
    return mFactory.make(contig, 0, mKmerSize);
  }


  List<Long> getLinks(boolean end, Kmer k, Map<Kmer, Long> contigEnds) {
    final List<Long> startLinks = new ArrayList<>();
    final long id;
    if (contigEnds.containsKey(k)) {
      id = contigEnds.get(k);
    } else if (contigEnds.containsKey(k.reverse())) {
      id = contigEnds.get(k.reverse());
    } else {
      return startLinks;
    }
    final Contig nextContig  = mContigGraph.contig(id);
    if (end) {
      if (startKmer(nextContig).equals(k)) {
        startLinks.add(id);
      }
      if (endKmer(nextContig).reverse().equals(k)) {
        startLinks.add(-id);
      }
    } else {
      if (startKmer(nextContig).reverse().equals(k)) {
        startLinks.add(-id);
      }
      if (endKmer(nextContig).equals(k)) {
        startLinks.add(id);
      }
    }
    return startLinks;
  }

  /**
   * Head along a contig in a the specified direction
   *
   *
   * @param sb Accumulate contig information into this
   * @param direction true means head right, false means head left
   * @param node starting Kmer node
   * @param visited set of nodes already visited while walking this contig
   * @throws IllegalStateException if we end up in a the middle of another contig
   */
  void walk(PreContig sb, boolean direction, Kmer node, Set<Kmer> visited) {
    Kmer current = node;
    while (true) {
      final Kmer next = uniqNext(current,  direction);
      if (next == null) {
        return;
      }
      if (uniqNext(next, !direction) == null) {
        return;
      }
      if (visited.contains(next)) {
        return;
      }
      visited.add(next);
      assert mDeBruijnGraph.contains(next);
      mDeBruijnGraph.setBuilt(next, true);
      current = next;
      sb.extend(direction, next, mDeBruijnGraph.frequency(next));

    }
  }

  /**
   * @param hash current Kmer
   * @param direction true means head right, false means head left
   * @return the unique Kmer that is one base in direction, or null if none exists or there are multiple
   */
  Kmer uniqNext(Kmer hash, boolean direction) {
    Kmer next = null;
    for (byte b = (byte) DNA.A.ordinal(); b <= DNA.T.ordinal(); ++b) {
      final Kmer successor = direction ? hash.successor(b) : hash.predecessor(b);
      if (next != null && mDeBruijnGraph.contains(successor)) {
        return null;
      } else if (mDeBruijnGraph.contains(successor)) {
        next = successor;
      }
    }
    return next;
  }
  static String minimalVersion(String key) {
    final String hash;
    final String complement = complement(key);
    if (key.compareTo(complement) < 0) {
      hash = key;
    } else {
      hash = complement;
    }
    return hash;
  }

  Pair<IntChunks, IntChunks> deleteTips(Pair<IntChunks, IntChunks> tipValues) {
    int tipCount = 0;
    long tipLength = 0;
    for (long i = 1; i <= mContigGraph.numberContigs(); ++i) {
      if (!mContigGraph.contigDeleted(i) && !tipValueOk(i, tipValues.getA(), tipValues.getB())) {
        tipLength += mContigGraph.contigLength(i);
        ++tipCount;
        mContigGraph.deleteContig(i);
      }
    }
    Diagnostic.info("Tip preContigs removed: " + tipCount);
    Diagnostic.info("Total length of tips: " + tipLength);
    return tipValues;
  }


  Pair<IntChunks, IntChunks> calculateTipValues() {
    //TODO Move me I don't belong in this class because I don't use the de Bruijn graph
    final IntChunks startTips = new IntChunks(mContigGraph.numberContigs() + 1);
    final IntChunks endTips = new IntChunks(mContigGraph.numberContigs() + 1);
    int changed = 1;
    while (changed > 0) {
      changed = 0;
      for (long i = 1; i <= mContigGraph.numberContigs(); ++i) {
        changed += updateTipValue(i, startTips, endTips) ? 1 : 0;
      }
    }
    return Pair.create(startTips, endTips);
  }


  boolean tipValueOk(long i, CommonIndex startTips, CommonIndex endTips) {
    final int tipThreshold = tipThreshold();
    final long startTip = startTips.get(i);
    if (startTip != 0 && startTip < tipThreshold) {
      return false;
    }
    final long endTip = endTips.get(i);
    return !(endTip != 0 && endTip < tipThreshold);
  }

  private int tipThreshold() {
    return mKmerSize * 2 + mTipConstant;
  }

  private boolean updateTipValue(long i, CommonIndex startTips, CommonIndex endTips) {
    int startTip = 0;
    int endTip = 0;
    final PathsIterator pathsIterator = mContigGraph.paths(i);
    long pathId;
    boolean startNull = false;
    boolean endNull = false;
    int startCount = 0;
    int endCount = 0;
    while ((pathId = pathsIterator.nextPathId()) != 0) {

      final int index = pathsIterator.contigIndex();
      final long linkedContig = mContigGraph.pathContig(pathId, 1 - index);
      if (index == 0) {
        ++endCount;
        final CommonIndex connected = linkedContig > 0 ? endTips : startTips;
        final long value = connected.get(Math.abs(linkedContig));
        if (value == 0) {
          endNull = true;
        } else {
          endTip = (int) Math.max(endTip, value);
        }
      } else {
        ++startCount;
        final CommonIndex connected = linkedContig > 0 ? startTips : endTips;
        final long value = connected.get(Math.abs(linkedContig));
        if (value == 0) {
          startNull = true;
        } else {
          startTip = (int) Math.max(startTip, value);
        }
      }
    }
    final int nextStartTip =  startCount == 0 ? mContigGraph.contigLength(i) : startTip + mContigGraph.contigLength(i) - mKmerSize + 1;
    final int nextEndTip = endCount == 0 ? mContigGraph.contigLength(i) : endTip + mContigGraph.contigLength(i) - mKmerSize + 1;
    boolean changed = false;
    if (!startNull && nextStartTip > startTips.get(i)) {
      changed = true;
      startTips.set(i, nextStartTip);
    }
    if (!endNull && nextEndTip > endTips.get(i)) {
      changed = true;
      endTips.set(i, nextEndTip);
    }
    return changed;
  }

}

