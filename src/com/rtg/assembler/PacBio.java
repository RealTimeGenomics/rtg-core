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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.UUID;

import com.rtg.alignment.ActionsHelper;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.ParamsTask;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.store.StoreDirProxy;

/**
 * Aligns PacBio reads data within a graph with the aim of resolving complicated sections.
 * We don't attempt to align things that look as if they will fall within a single contig
 */
public final class PacBio extends ParamsTask<PacBioParams, PacBioStatistics> {

  static final double ERROR_RATE = 0.06;
  static final int ERROR_FLOOR = 3;

  private static final int WORDSIZE = 15;
  private static final int LINE_LENGTH = 120;

  /**
   * Construct a pac bio aligner
   * @param params Params object
   * @param out somewhere to dump output
   */
  public PacBio(PacBioParams params, OutputStream out) {
    super(params, out, new PacBioStatistics(params.directory()), null);
  }

  void mapPacBio(ReadPairSource source, MutableGraph graph) throws IOException {
    final GraphIndex index = new GraphIndex(graph, 1, WORDSIZE);
    List<byte[]> fragments;
    final PathTracker tracker = new PathTracker(new PalindromeTracker(graph));
    while ((fragments = source.nextFragments()) != null) {
      for (byte[] read : fragments) {
        final List<List<ContigPosition>> hits = index.hits(read, graph, index.getSearchFunction());
        boolean hasHit = false;
        for (List<ContigPosition> hit : hits) {
          //out.println("Dodgy read: " + readId);
          //printHits(hits);
          if (hit.size() > 0) {
            hasHit = true;
            break;
          }
        }
        if (!hasHit) {
          continue;
        }
        final boolean isInternal = isInternal(graph, hits);
        if (isInternal) {
          mStatistics.increment(PacBioStatistics.Stat.INTERNAL_READS);
        }
        if (!isInternal) {
          //printHits(hits);
          final HitMap mergedHits = joinHits(hits, WORDSIZE);
          if (mergedHits.size() > 1) {
            final List<PartialAlignment> alignments = alignHits(mergedHits, graph, read);
            if (alignments.size() >= 2) {
              Collections.sort(alignments);
              int endPoint = 0;
              for (PartialAlignment pa : alignments) {
                endPoint = Math.max(pa.getReadEnd(), endPoint);
              }
              final SortedSet<PartialAlignment> sorted = new TreeSet<>();
              sorted.addAll(alignments);
              final Map<Long, List<PacBioPath>> longListMap = joinAlignments(new ArrayList<>(sorted), graph);
              final PacBioPath best = uniqueBest(longListMap, read.length);
              if (best != null) {

                mStatistics.increment(PacBioStatistics.Stat.CROSS_CONTIG);
                tracker.increment(best.toPath());
              }
            }
          }
        }
      }
      mStatistics.increment(PacBioStatistics.Stat.TOTAL_READS);
    }
    final SortedMap<List<Long>, Integer> merge = PathTracker.merge(Collections.singletonList(tracker));
    PathTracker.apply(merge, graph);
  }

  static boolean isInternal(Graph graph, List<List<ContigPosition>> hits) {
    final int[] hitCounts = new int[(int) graph.numberContigs() + 1];
    final int hitsSize = hits.size();
    for (final List<ContigPosition> positionHits : hits) {
      for (ContigPosition pos : positionHits) {
        hitCounts[(int) Math.abs(pos.mContigId)]++;
      }
    }
    int maxIndex = 0;
    for (int i = 0; i < hitCounts.length; i++) {
      if (hitCounts[i] > hitCounts[maxIndex]) {
        maxIndex = i;
      }
    }
    if (maxIndex == 0) {
      return false;
    }

    final int contigLength = graph.contigLength(maxIndex);
    boolean isInternal = false;
    if (contigLength >= hits.size()) {
      isInternal = true;
      for (int i = 0; i < hitsSize; i++) {
        final List<ContigPosition> positionHits = hits.get(i);
        for (ContigPosition pos : positionHits) {
          if (Math.abs(pos.mContigId) == maxIndex) {
            if (i > pos.mPosition) {
              isInternal = false;
            }
            if (hits.size() - i > contigLength - pos.mPosition) {
              isInternal = false;
            }
          }
        }
      }
    }
    return isInternal;
  }

  static Map<Long, List<PacBioPath>> joinAlignments(List<PartialAlignment> alignments, MutableGraph graph) {
    final Map<Long, List<PacBioPath>> pathMap = new HashMap<>();
    for (PartialAlignment part : alignments) {
      final Set<Long> predecessors = MergeNodes.predecessors(graph, part.getContig());
      boolean foundPrevious = false;
      final List<PacBioPath> newPathList = new ArrayList<>();
      for (long previousContig : predecessors) {
        final List<PacBioPath> paths = pathMap.get(previousContig);
        if (paths != null) {
          for (PacBioPath path : paths) {
            if (Math.abs(path.mAlignment.getReadEnd() - part.getReadStart() - graph.contigOverlap()) < graph.contigOverlap() / 2) {
              newPathList.add(new PacBioPath(path, part));
              foundPrevious = true;
              path.mIsPrefix = true;
            }
          }
        }
      }
      final List<PacBioPath> mapPaths = getOrAdd(pathMap, part.getContig());
      for (PacBioPath newPath : newPathList) {
//        out.println("considering: " + newPath.toString());
        addOrReplace(mapPaths, newPath);
      }
      if (!foundPrevious && part.getReadStart() == 0) {
        final List<PacBioPath> addingPaths = getOrAdd(pathMap, part.getContig());
        addingPaths.add(new PacBioPath(null, part));
      }
    }
    return pathMap;
  }

  private static void addOrReplace(List<PacBioPath> mapPaths, PacBioPath pacBioPath) {
    final List<PacBioPath> remove = new ArrayList<>();
    boolean better = true;
    final int newScore = pacBioPath.score();
    final long startContig = pacBioPath.mFirstContig;
    for (PacBioPath existing : mapPaths) {
      final int existingScore = existing.score();
      final long existingStartContig = existing.mFirstContig;
      if (existingStartContig == startContig && existing.mReadStart == pacBioPath.mReadStart && pacBioPath.mAlignment.getReadEnd() == existing.mAlignment.getReadEnd()) {
        // match found
        if (newScore < existingScore) {
          remove.add(existing);
        } else if (newScore == existingScore) {
          if (!existing.toPath().equals(pacBioPath.toPath())) {
            existing.mIsDupe = true;
          }
          better = false;
        } else {
          better = false;
        }

      }

    }
    mapPaths.removeAll(remove);
    if (better) {
      mapPaths.add(pacBioPath);
    }
  }

  static PacBioPath uniqueBest(Map<Long, List<PacBioPath>> pathMap, int readLength) {
    final List<PacBioPath> allPaths = new ArrayList<>();
    for (List<PacBioPath> pathList : pathMap.values()) {
      allPaths.addAll(pathList);
    }
    final List<PacBioPath> wholeReadPaths = new ArrayList<>();
    for (PacBioPath path : allPaths) {
      if (path.mIsPrefix || path.mReadStart > 0 || path.mAlignment.getReadEnd() != readLength) {
        continue;
      }
      wholeReadPaths.add(path);
    }
    Collections.sort(wholeReadPaths, new ScoreComparator());
    if (wholeReadPaths.size() > 0) {
      final PacBioPath bestPath = wholeReadPaths.get(0);
      if ((wholeReadPaths.size() == 1 || wholeReadPaths.get(1).score() > bestPath.score()) && !bestPath.mIsDupe) {
        if (bestPath.mPrevious != null) {
          return bestPath;
        }
      }
    }
    return null;
  }

  @Override
  protected void exec() throws IOException {
    final MutableGraph g = (MutableGraph) GraphReader.read(new StoreDirProxy(mParams.graph()));
    if (mParams.trimGraph()) {
      final int trimmed = GraphCleanup.clean(GraphCleanup.MIN_LENGTH, g);
      Diagnostic.info("Trimmed " + trimmed + " short contigs");
    }
    for (File readDir : mParams.reads()) {
      try (final ReadPairSource source = new ReadPairSource(SequencesReaderFactory.createDefaultSequencesReader(readDir))) {
        mapPacBio(source, g);
      }
    }
    final Graph sortedGraph = GraphSorter.sortedGraph(g);
    GraphWriter.write(sortedGraph, new StoreDirProxy(mParams.directory()), "pac_bio_map", Collections.<UUID>emptySet());
  }

  static class ScoreComparator implements Comparator<PacBioPath>, Serializable {
    @Override
    public int compare(PacBioPath o1, PacBioPath o2) {
      return Integer.valueOf(o1.score()).compareTo(o2.score());
    }
  }

  private static List<PacBioPath> getOrAdd(Map<Long, List<PacBioPath>> pathMap, long contig) {
    List<PacBioPath> newPathList = pathMap.get(contig);
    if (newPathList == null) {
      newPathList = new ArrayList<>();
      pathMap.put(contig, newPathList);
    }
    return newPathList;
  }

  static List<PartialAlignment> alignHits(final HitMap hits, final MutableGraph graph, final byte[] read) {
    final List<PartialAlignment> alignments = new ArrayList<>();
    for (Map.Entry<Long, List<HitCollection>> longListEntry : hits.entrySet()) {
      final List<HitCollection> value = longListEntry.getValue();
      for (final HitCollection next : value) {
        final PartialAlignment alignForward = alignHit(next, graph, read);
        if (alignForward != null) {
          final PartialAlignment alignBackward = alignHitRightFixed(next, graph, read);
          if (alignBackward != null) {
            final PartialAlignment joined = new PartialAlignment(alignForward.getAlignmentScore() + alignBackward.getAlignmentScore(), alignBackward.getReadStart(), alignForward.getReadEnd(), alignForward.getContig(), alignBackward.getContigStart(), alignForward.getContigEnd());
            alignments.add(joined);
          }
        }
      }
    }
    return alignments;
  }

  private static int[] actionsCopy(final int[] actions) {
    return Arrays.copyOf(actions, actions.length);
  }

  static PartialAlignment alignHit(final HitCollection hc, final MutableGraph graph, final byte[] read) {
    final long contigId = hc.mContigId;
    final HitPosition hitPosition = hc.mPositions.get(0);
    final int contigPos = hitPosition.contigPosition();
    final byte[] contig = getContigArray(graph, contigId);
    final PathAligner pa = new PathAligner(graph);
    final int readPos = hitPosition.readPosition();
    final int shift = WORDSIZE - 1;
    final int rStart = readPos - shift;
    final int cStart = contigPos - shift;

    final int length = Math.min(read.length - rStart, contig.length - cStart);
    PartialAlignment partial = null;
    if (read.length - rStart <= contig.length - cStart) {

      final int[] actionsReadFixed = actionsCopy(pa.alignRight(read, rStart, read.length, contig, cStart));
      final int alignmentScoreReadFixed = ActionsHelper.alignmentScore(actionsReadFixed);
//      System.err.println("one");
      if (alignmentScoreReadFixed < length) {
        final int alignmentLengthContig = ActionsHelper.actionsCount(actionsReadFixed) - ActionsHelper.deletionFromReadAndOverlapCount(actionsReadFixed);
        final int alignmentLengthRead = ActionsHelper.actionsCount(actionsReadFixed) - ActionsHelper.insertionIntoReadAndGapCount(actionsReadFixed);
        partial = new PartialAlignment(alignmentScoreReadFixed, rStart, rStart + alignmentLengthRead, contigId, cStart, cStart + alignmentLengthContig);
      }
    } else {
//      out.println("Contig: " + DnaUtils.bytesToSequenceIncCG(contig));
//      out.println("Read: " + DnaUtils.bytesToSequenceIncCG(read));
//      out.println(hitPosition);
      final int[] actionsContigFixed = actionsCopy(pa.alignRight(contig, cStart, contig.length, read, rStart));
      final int alignmentScoreContigFixed = ActionsHelper.alignmentScore(actionsContigFixed);
//      System.err.println("two");
      if (alignmentScoreContigFixed < length) { // && alignmentScoreContigFixed < alignmentScoreReadFixed) {
        final int alignmentLengthRead = ActionsHelper.actionsCount(actionsContigFixed) - ActionsHelper.deletionFromReadAndOverlapCount(actionsContigFixed);
        final int alignmentLengthContig = ActionsHelper.actionsCount(actionsContigFixed) - ActionsHelper.insertionIntoReadAndGapCount(actionsContigFixed);
        partial = new PartialAlignment(alignmentScoreContigFixed, rStart, rStart + alignmentLengthRead, contigId, cStart, cStart + alignmentLengthContig);
      }
    }

    return partial;
  }

  private static byte[] getContigArray(MutableGraph graph, long contigId) {
    final byte[] contigNt = new byte[graph.contigLength(contigId)];
    for (int k = 0; k < contigNt.length; k++) {
      contigNt[k] = graph.nt(contigId, k);
    }
    return contigNt;
  }

  static PartialAlignment alignHitRightFixed(final HitCollection hc, final MutableGraph graph, final byte[] read) {
    final long contigId = hc.mContigId;
    final HitPosition hitPosition = hc.mPositions.get(0);
    final int contigPos = hitPosition.contigPosition();
    final byte[] contigNt = getContigArray(graph, contigId);
    final PathAligner pa = new PathAligner(graph);
    final int readPos = hitPosition.readPosition();

    // Treat the longer of read, contig as the template in the aligner.  The aligner can tolerate
    // some shift due to running off the end, but data structure size during alignment is
    // expensive if read size is too large.
    PartialAlignment partial = null;
    if (readPos <= contigPos) {
      final int[] actionsReadFixed = pa.alignLeft(read, 0, readPos, contigNt, contigPos);
      final int alignmentScoreReadFixed = ActionsHelper.alignmentScore(actionsReadFixed);
      final int contigStart = ActionsHelper.zeroBasedTemplateStart(actionsReadFixed);
      final int readStart = readPos - ActionsHelper.actionsCount(actionsReadFixed) + ActionsHelper.insertionIntoReadAndGapCount(actionsReadFixed);
//      System.err.println(contigPos + " " + ActionsHelper.actionsCount(actionsReadFixed) + " " + ActionsHelper.deletionFromReadCount(actionsReadFixed));
//      System.err.println("three");
      // Score should surely be less than the portion aligned
      if (alignmentScoreReadFixed < readPos) {
        partial = new PartialAlignment(alignmentScoreReadFixed, readStart, readPos, contigId, contigStart, contigPos);
      }
    } else {
      final int[] actionsContigFixed = pa.alignLeft(contigNt, 0, contigPos, read, readPos);
      final int alignmentScoreContigFixed = ActionsHelper.alignmentScore(actionsContigFixed);
      final int readStart = readPos - ActionsHelper.actionsCount(actionsContigFixed) + ActionsHelper.deletionFromReadAndOverlapCount(actionsContigFixed);
      final int contigStart = contigPos - ActionsHelper.actionsCount(actionsContigFixed) + ActionsHelper.insertionIntoReadAndGapCount(actionsContigFixed);
//      System.err.println("four");
      if (alignmentScoreContigFixed < contigPos) {
        partial = new PartialAlignment(alignmentScoreContigFixed, readStart, readPos, contigId, contigStart, contigPos);
      }
    }
    return partial;
  }

  static void printAlignment(final String header, final byte[] read, final byte[] contigNt, final int rStart, final int cStart, final String a, PrintStream out) {
    final StringBuilder rOut = new StringBuilder();
    final StringBuilder tOut = new StringBuilder();
    for (int k = 0, r = rStart, t = cStart; k < a.length(); k++) {
      switch (a.charAt(k)) {
        case 'X':
        case '=':
          rOut.append(DnaUtils.getBase(read[r++]));
          tOut.append(DnaUtils.getBase(contigNt[t++]));
          break;
       case 'D':
          rOut.append(DnaUtils.getBase(read[r++]));
          tOut.append("-");
          break;
       case 'I':
          rOut.append("-");
          tOut.append(DnaUtils.getBase(contigNt[t++]));
          break;
       default:
            throw new RuntimeException();
      }
    }
    out.println(header);
    final String tA = tOut.toString();
    final String rA = rOut.toString();
    for (int k = 0; k < a.length(); k += LINE_LENGTH) {
      out.println(tA.substring(k, Math.min(k + LINE_LENGTH, a.length())));
      out.println(a.substring(k, Math.min(k + LINE_LENGTH, a.length())));
      out.println(rA.substring(k, Math.min(k + LINE_LENGTH, a.length())));
      out.println();
    }
  }


  static class HitMap extends HashMap<Long, List<HitCollection>> {
    List<HitCollection> getOrAdd(Long name) {
      if (containsKey(name)) {
        return get(name);
      }
      final List<HitCollection> res = make();
      put(name, res);
      return res;
    }
    List<HitCollection> make() {
      return new ArrayList<>();
    }

  }

  static HitMap joinHits(List<List<ContigPosition>> hits, int wordSize) {
    final HitMap map = new HitMap();
    final int size = hits.size();
    for (int readPosition = 0; readPosition < size; readPosition++) {
      final List<ContigPosition> posList = hits.get(readPosition);
      for (ContigPosition pos : posList) {
        final List<HitCollection> list = map.getOrAdd(pos.mContigId);
        int bestMatch = -1;
        int bestOffset = Integer.MAX_VALUE;
        int bestDistance = 0;
        boolean ambiguous = false;

        final int listSize = list.size();
        for (int i = 0; i < listSize; i++) {
          final HitCollection collection = list.get(i);
          final int readDistance = readPosition - collection.top().readPosition();
          final int contigDistance = pos.mPosition - collection.top().contigPosition();
          final int offset = Math.abs(contigDistance - readDistance);
          if (offset != 0 && (contigDistance < wordSize || readDistance < wordSize)) {
            continue;
          }
          if (offset < bestOffset) {
            bestMatch = i;
            bestOffset = offset;
            ambiguous = false;
            bestDistance =  contigDistance;
          } else if (offset == bestOffset) {
            ambiguous = true;
          }
        }
        if (bestOffset < Math.max(bestDistance * ERROR_RATE, ERROR_FLOOR) && !ambiguous) {
          list.get(bestMatch).add(new HitPosition(pos.mPosition, readPosition));
        } else {
          final HitCollection collect = new HitCollection(pos.mContigId);
          collect.add(new HitPosition(pos.mPosition, readPosition));
          list.add(collect);
        }
      }
    }
    return map;
  }
  static class HitCollection {
    long mContigId;
    List<HitPosition> mPositions = new LinkedList<>();
    HitCollection(long contigId) {
      mContigId = contigId;
    }
    HitPosition top() {
      return mPositions.get(mPositions.size() - 1);
    }
    void add(HitPosition pos) {
      mPositions.add(pos);
    }

    @Override
    public String toString() {
      return "HitCollection" + mPositions;
    }
  }

  static class HitPosition {
    final int mContigPosition;
    final int mReadPosition;
    HitPosition(int contig, int read) {
      mContigPosition = contig;
      mReadPosition = read;
    }
    int readPosition() {
      return mReadPosition;
    }
    int contigPosition() {
      return mContigPosition;
    }
    public String toString() {
      return "Hit{C" + mContigPosition + ", R" + mReadPosition + "}";
    }

  }
}
