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
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.Histogram;
import com.rtg.util.IORunnable;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MathUtils;
import com.rtg.util.NormalDistribution;
import com.rtg.util.ProgramState;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.Utils;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class GraphMap {

  static final int MAX_HITS_PER_START_POSITION = GlobalFlags.getIntegerValue(CoreGlobalFlags.ASSEMBLER_MAX_HITS_PER_START_POS_FLAG);
  static final int INSERT_DEVIATIONS = GlobalFlags.getIntegerValue(CoreGlobalFlags.ASSEMBLER_INSERT_DEVIATIONS_FLAG);

  private final GraphIndex mIndex;

  PairConstraintWriter mOut;
  private final GraphMapStatistics mStatistics;
  private final MutableGraph mGraph;
  private final IntChunks mContigReadCounts;
  private final IntChunks mPathReadCounts;

  public PathTracker getPathTracker() {
    return mPathTracker;
  }

  private final PathTracker mPathTracker;

  public GraphMapStatistics getStatistics() {
    return mStatistics;
  }
  public ConstraintCache getConstraints() {
    return mOut.mCache;
  }

  GraphMap(GraphIndex index, MutableGraph graph, PairConstraintWriter output, PathTracker pathTracker) {
    mIndex = index;
    mGraph = graph;
    mStatistics = new GraphMapStatistics(null);
    mContigReadCounts = new IntChunks(graph.numberContigs() + 1);
    mPathReadCounts = new IntChunks(graph.numberPaths() + 1);
    mPathTracker = pathTracker;
    mOut = output;
  }

  static long findPath(List<Long> path, Graph graph, boolean palindromesEqual) {
    final long firstContig = path.get(0);
    final boolean firstIsPalindrome = GraphAlignment.isPalindrome(firstContig, graph);
    final long first = firstIsPalindrome ? Math.abs(firstContig) : firstContig;
    final PathsIterator iterator = graph.paths(first);
    final long result = findPathInternal(path, iterator, graph, palindromesEqual);
    if (result == 0 && palindromesEqual && firstIsPalindrome) {
      final PathsIterator reverseIterator = graph.paths(-first);
      return findPathInternal(path, reverseIterator, graph,  palindromesEqual);
    }
    return result;
  }

  static long findPathInternal(List<Long> path, PathsIterator iterator, Graph graph, boolean palindromesEqual) {
    long pathId;
    pathsLoop: while ((pathId = iterator.nextPathId()) != 0) {
      if (graph.pathLength(pathId) != path.size()) {
        continue;
      }
      for (int i = 0; i < path.size(); ++i) {
        final Long contig = path.get(i);
        if (contig != graph.pathContig(pathId, i) && !(palindromesEqual && contig == -graph.pathContig(pathId, i) && GraphAlignment.isPalindrome(contig, graph))) {
          continue pathsLoop;
        }
      }
      return pathId;
    }
    return 0;
  }

  static void setInsertSizes(List<ReadPairSource> reads, GraphMapParams params, MutableGraph graph, GraphIndex index) throws IOException {
    final SimpleThreadPool pool = new SimpleThreadPool(params.numberThreads(), "InsertCalc", true);
    final List<ReadPairSource> paired = new ArrayList<>();
    for (final ReadPairSource reader : reads) {
      reader.reset();
      if (reader.numberFragments() > 1) {
        paired.add(reader);
      }
    }
    final AsyncReadPool readPool = new AsyncReadPool("ReadForInsertCalclation", paired);
    final List<InsertSizeRunnable> insertJobs = new ArrayList<>();
    for (int i = 0; i < params.numberThreads(); ++i) {
      final GraphMap graphMap = new GraphMap(index, graph, null, new PathTracker(new PalindromeTracker(graph)));
      final InsertSizeRunnable run = new InsertSizeRunnable(graphMap, params.maxMismatches(), readPool.sources());
      insertJobs.add(run);
      pool.execute(run);
    }
    pool.terminate();
    readPool.terminate();
    for (final ReadPairSource reader : paired) {
      final List<NormalDistribution> deviations = new ArrayList<>();
      for (InsertSizeRunnable job : insertJobs) {
        final NormalDistribution e = job.mDeviationMap.get(reader);
        deviations.add(e);
        Diagnostic.developerLog(e.toString());
      }
      final NormalDistribution combined = new NormalDistribution();
      combined.add(deviations);
      Diagnostic.info("Insert size distribution: " + combined);
      reader.setMinInsertSize((int) MathUtils.round(combined.mean() - combined.stdDev() * INSERT_DEVIATIONS));
      Diagnostic.info("Min Insert: " + reader.minInsertSize());
      reader.setMaxInsertSize((int) MathUtils.round(combined.mean() + combined.stdDev() * INSERT_DEVIATIONS));
      Diagnostic.info("Max Insert: " + reader.maxInsertSize());
    }
  }

  private void progress(final long id, final long total) {
    if (id % 100000 == 0 || id == total) {
      Diagnostic.progress(id + "/" + total + " (" + Utils.realFormat(100.0 * id / total, 1) + "%)");
    }
  }

  void mapReads(AsyncReadSource reader, IntegerOrPercentage mismatches) throws IOException {
    final GraphTraversions traverse = new GraphTraversions(mGraph);
    final GraphAligner aligner = reader.aligner(mGraph, mismatches, traverse);
    final PairJoiner joiner = new PairJoiner(mGraph, getOverlap(), traverse);
    final Histogram asHistogram = new Histogram();
    final long numberReads = reader.getNumberReads();
    final AlignmentIterator it = new AlignmentIterator(reader, mGraph, aligner, mIndex, mStatistics);
    final int expectedInsert = (reader.minInsertSize() + reader.maxInsertSize()) / 2;
    while (it.hasNext()) {
      ProgramState.checkAbort();
      final AlignmentIterator.ReadAlignment alignment = it.next();
      progress(alignment.mId, numberReads);
      final List<Set<GraphAlignment>> fragmentAlignments = alignment.mFragments;
      if (fragmentAlignments.size() > 1) {
        mStatistics.increment(GraphMapStatistics.Stat.PAIRED_END);
        final Set<GraphAlignment> pairAlignments = joiner.paired(fragmentAlignments, reader.minInsertSize(), reader.maxInsertSize());
        if (pairAlignments == PairJoiner.TOO_MANY_PATHS) {
          mStatistics.increment(GraphMapStatistics.Stat.TOO_MANY_PAIR_PATHS);
        } else if (pairAlignments.size() == 0) {
          mStatistics.increment(GraphMapStatistics.Stat.NO_PAIRINGS);
        } else if (addBestAlignment(pairAlignments, true, asHistogram) != null) {
          mStatistics.increment(GraphMapStatistics.Stat.PAIRED);
          continue;
        }
      } else {
        mStatistics.increment(GraphMapStatistics.Stat.SINGLE_END);
      }
      boolean singleBest = true;

      final List<GraphAlignment> bestList = new ArrayList<>();
      for (Set<GraphAlignment> alignments : fragmentAlignments) {
        final GraphAlignment best = addBestAlignment(alignments, false, asHistogram);
        singleBest &= best != null;
        bestList.add(best);
      }
      if (singleBest && bestList.size() == 2) {
        final GraphAlignment first = bestList.get(0);
        final GraphAlignment second = bestList.get(1);
        final long contig1 = first.endContig();
        final long contig2 = second.endContig();
        final boolean leftEnd = traverse.get(contig1).mNext.size() < 1;
        final boolean rightEnd = traverse.get(contig2).mNext.size() < 1;

        final int distanceFromEnd1 = mGraph.contigLength(contig1) - first.endPosition() - 1;
        final int distanceFromEnd2 = mGraph.contigLength(contig2) - second.endPosition() - 1;
        if (leftEnd && rightEnd && distanceFromEnd1 + distanceFromEnd2 < reader.maxInsertSize()) {
          //System.err.println("contig1:" + contig1 + " contig2:" +  contig2 + " distance1:" + distanceFromEnd1 + "distance2:" + distanceFromEnd2 + "expectedInsert: " + expectedInsert);
          mOut.mCache.addConstraint(contig1, contig2, distanceFromEnd1, distanceFromEnd2, expectedInsert);
        }
        mOut.writeConstraint(contig1, contig2, distanceFromEnd1, distanceFromEnd2, leftEnd, rightEnd);
      }
    }
    progress(numberReads, numberReads);
  }

  private GraphAlignment addBestAlignment(Set<GraphAlignment> alignments, boolean paired, Histogram asHistogram) {
    GraphAlignment best = null;
    int score = Integer.MAX_VALUE;
    int scoreCount = 0;
    for (GraphAlignment alignment : alignments) {
      if (alignment.mScore < score) {
        best = alignment;
        score = alignment.mScore;
        scoreCount = 1;
      } else if (alignment.mScore == score) {
        ++scoreCount;
      }
    }
    if (scoreCount > 1) {
      mStatistics.increment(paired ? GraphMapStatistics.Stat.TOO_MANY_PAIR_PATHS : GraphMapStatistics.Stat.TOO_MANY_PATHS);
      return null;
    } else if (best == null) {
      mStatistics.increment(GraphMapStatistics.Stat.NO_PATHS);
      return null;

    }
    asHistogram.increment(score);
    if (!paired) {
      mStatistics.increment(GraphMapStatistics.Stat.MAPPED);
    }
    if (best.contigs().size() == 1) {
//      mOut.println(readId + "\tC" + best.contigs().get(0) + "\t" + best.startPosition() + "\t" + best.endPosition());
      mStatistics.increment(paired ? GraphMapStatistics.Stat.SINGLE_CONTIG_PAIR : GraphMapStatistics.Stat.SINGLE_CONTIG_MAPPING);
      singleContig(best);
    } else {
      mStatistics.increment(paired ? GraphMapStatistics.Stat.CROSS_CONTIG_PAIR : GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE);
      crossContig(best);
    }
    return best;
  }

  private void crossContig(GraphAlignment best) {
    mPathTracker.increment(best.contigs());
  }

  static int existingReadCount(MutableGraph graph, long pathId) {
    final String readCount = graph.pathAttribute(pathId, GraphKmerAttribute.READ_COUNT);
    final int oldCount;
    if (readCount == null) {
      oldCount = 0;
    } else {
      oldCount = Integer.parseInt(readCount);
    }
    return oldCount;
  }

  private void singleContig(GraphAlignment best) {
    final long contigId = best.contigs().get(0);
    final long absId = Math.abs(contigId);
    mContigReadCounts.set(absId, mContigReadCounts.get(absId) + 1);
  }

  /**
   * Calculate the insert size of a set of reads from reads that map to the same contig
   * @param reader source of paired reads
   * @param mismatches the number of mismatches allowed in alignment
   * @throws IOException when reading the input data
   * @return the distribution calculated for the insert sizes
   */
  public NormalDistribution calculateInserts(AsyncReadSource reader, IntegerOrPercentage mismatches) throws IOException {
    assert reader.getNumberFragments() > 1;
    final GraphAligner aligner = reader.aligner(mGraph, mismatches, new NullTraversions());
    final AlignmentIterator it = new AlignmentIterator(reader, mGraph, aligner, mIndex, new GraphMapStatistics(null));
    final NormalDistribution stdDev = new NormalDistribution();
    while (it.hasNext()) {
      final AlignmentIterator.ReadAlignment alignment = it.next();
      final List<Set<GraphAlignment>> fragmentAlignments = alignment.mFragments;
      if (fragmentAlignments.size() > 1) {
        final Set<GraphAlignment> leftSet = fragmentAlignments.get(0);
        final Set<GraphAlignment> rightSet = fragmentAlignments.get(1);
        if (leftSet.size() == 1 && rightSet.size() == 1) {
          final GraphAlignment leftAlignment = setToVal(leftSet);
          final GraphAlignment rightAlignment = setToVal(rightSet);
          if (leftAlignment.contigs().size() == 1 && rightAlignment.contigs().size() == 1 && leftAlignment.startContig() == -rightAlignment.startContig()) {
            final long insert = mGraph.contigLength(rightAlignment.startContig()) - rightAlignment.endPosition() - 1 - leftAlignment.endPosition();
            stdDev.add(insert);
          }
        }
      }
    }
    return stdDev;
  }

  <T> T setToVal(Set<T> set) {
    for (T val : set) {
      return val;
    }
    return null;
  }

  static void finalizeCounts(List<GraphMap> mappings, MutableGraph graph) {
    for (long i = 1; i < graph.numberContigs() + 1; ++i) {
      final String readCount = graph.contigAttribute(i, GraphKmerAttribute.READ_COUNT);
      int count;
      if (readCount == null) {
        count = 0;
      } else {
        count = Integer.parseInt(readCount);
      }
      for (GraphMap mapping : mappings) {
        if (mapping.mContigReadCounts.length() > i) {
          count += mapping.mContigReadCounts.get(i);
        }
      }
      if (count > 0) {
        graph.setContigAttribute(i, GraphKmerAttribute.READ_COUNT, "" + count);
      }
    }
    for (long i = 1; i < graph.numberPaths() + 1; ++i) {
      final String readCount = graph.pathAttribute(i, GraphKmerAttribute.READ_COUNT);
      int count;
      if (readCount == null) {
        count = 0;
      } else {
        count = Integer.parseInt(readCount);
      }
      for (GraphMap mapping : mappings) {
        if (mapping.mPathReadCounts.length() > i) {
          count += mapping.mPathReadCounts.get(i);
        }
      }
      if (count > 0) {
        graph.setPathAttribute(i, GraphKmerAttribute.READ_COUNT, "" + count);
      }
    }
  }

  public int getOverlap() {
    return mGraph.contigOverlap();
  }

  static class MapRunnable implements IORunnable {
    final GraphMap mMapper;
    final IntegerOrPercentage mMismatches;
    final List<AsyncReadSource> mReaders;

    MapRunnable(GraphMap mapper, IntegerOrPercentage mismatches, List<AsyncReadSource> reader) {
      mMapper = mapper;
      mMismatches = mismatches;
      mReaders = reader;
    }

    @Override
    public void run() throws IOException {
      for (AsyncReadSource reader : mReaders) {
        mMapper.mapReads(reader, mMismatches);
      }
    }

  }
  static class InsertSizeRunnable implements IORunnable {
    final GraphMap mMapper;
    final IntegerOrPercentage mMismatches;
    final List<AsyncReadSource> mReaders;
    final Map<ReadPairSource, NormalDistribution> mDeviationMap = new HashMap<>();

    InsertSizeRunnable(GraphMap mapper, IntegerOrPercentage mismatches, List<AsyncReadSource> reader) {
      mMapper = mapper;
      mMismatches = mismatches;
      mReaders = reader;
    }

    @Override
    public void run() throws IOException {
      int readerId = 0;
      for (AsyncReadSource reader : mReaders) {
        if (reader.getNumberFragments() > 1) {
          mDeviationMap.put(reader.underlying(), mMapper.calculateInserts(reader, mMismatches));
        } else {
          mDeviationMap.put(reader.underlying(), new NormalDistribution());
        }
        Diagnostic.developerLog("Finished insertSize run for reader " + readerId++);
      }
    }
  }
}
