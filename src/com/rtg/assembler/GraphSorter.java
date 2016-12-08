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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.Path;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.util.store.StoreDirProxy;

/**
 */
final class GraphSorter {
  private GraphSorter() { }

  static class OriginalContig implements Comparable<OriginalContig> {
    OriginalContig(Graph g, long contigId) {
      final int limit = g.contigLength(contigId);
      final Contig forward = g.contig(contigId);
      final Contig reverse = g.contig(-contigId);
      boolean forwardIsSmaller = true;
      for (int i = 0; i < limit; ++i) {
        if (forward.nt(i) > reverse.nt(i)) {
          forwardIsSmaller = false;
          break;
        } else if (forward.nt(i) < reverse.nt(i)) {
          break;
        }
      }
      if (forwardIsSmaller) {
        mContig = forward;
        mOriginalId = contigId;
      } else {
        mContig = reverse;
        mOriginalId = -contigId;
      }
    }
    long mOriginalId;
    Contig mContig;

    @Override
    public int compareTo(OriginalContig o) {
      final int limit = Math.min(o.mContig.length(), mContig.length());
      for (int i = 0; i < limit; ++i) {
        final int compare = Byte.compare(mContig.nt(i), o.mContig.nt(i));
        if (compare != 0) {
          return compare;
        }
      }
      return Integer.compare(mContig.length(), o.mContig.length());
    }
    @Override
    public boolean equals(Object o) {
      return o != null && o.getClass() == this.getClass() && compareTo((OriginalContig) o) == 0;
    }

    @Override
    public int hashCode() {
      int result = (int) (mOriginalId ^ (mOriginalId >>> 32));
      result = 31 * result + (mContig != null ? mContig.hashCode() : 0);
      return result;
    }

    @Override
    public String toString() {
      return "OriginalContig{" + "mOriginalId=" + mOriginalId + '}';
    }
  }
  static class TranslatedPath implements Comparable<TranslatedPath> {
    TranslatedPath(Graph g, long pathId, NegativeMap translate) {
      final int limit = g.pathLength(pathId);
      final Path original = g.path(pathId);

      final long[] translated = new long[original.length()];
      final long[] reverseTranslated = new long[original.length()];
      for (int i = 0; i < limit; ++i) {
        final long contig = translate.get(original.contig(i));
        translated[i] = contig;
        reverseTranslated[limit - i - 1] = -contig;
      }
      final Path forward = new PathArray(translated);
      final Path reverse = new PathArray(reverseTranslated);

      boolean forwardIsSmaller = true;
      for (int i = 0; i < limit; ++i) {
        if (forward.contig(i) > reverse.contig(i)) {
          forwardIsSmaller = false;
          break;
        } else if (forward.contig(i) < reverse.contig(i)) {
          break;
        }
      }
      if (forwardIsSmaller) {
        mOriginalId = pathId;
        mPath = forward;
      } else {
        mOriginalId = -pathId;
        mPath = reverse;
      }
    }
    long mOriginalId;
    Path mPath;

    @Override
    public int compareTo(TranslatedPath o) {
      final int limit = Math.min(o.mPath.length(), mPath.length());
      for (int i = 0; i < limit; ++i) {
        final int compare = Long.compare(mPath.contig(i), o.mPath.contig(i));
        if (compare != 0) {
          return compare;
        }
      }
      return Integer.compare(mPath.length(), o.mPath.length());
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (o == null || getClass() != o.getClass()) {
        return false;
      }

      final TranslatedPath that = (TranslatedPath) o;
      return this.compareTo(that) == 0;

    }

    @Override
    public int hashCode() {
      int result = (int) (mOriginalId ^ (mOriginalId >>> 32));
      result = 31 * result + (mPath != null ? mPath.hashCode() : 0);
      return result;
    }
  }
  static class NegativeMap {
    Map<Long, Long> mMap = new HashMap<>();
    void put(long key, long value) {
      if (key < 0) {
        mMap.put(-key, -value);
      } else {
        mMap.put(key, value);
      }
    }
    long get(long key) {
      if (key < 0) {
        return -mMap.get(-key);
      } else {
        return mMap.get(key);
      }
    }
    @Override
    public String toString() {
      return mMap.toString();
    }
  }
  static Graph sortedGraph(Graph unsorted) {
    final GraphKmerAttribute sorted = new GraphKmerAttribute(unsorted.contigOverlap(), unsorted.contigAttributes(), unsorted.pathAttributes());
    final List<OriginalContig> sortMe = new ArrayList<>();
    for (long i = 1; i <= unsorted.numberContigs(); ++i) {
      if (!unsorted.contigDeleted(i)) {
        sortMe.add(new OriginalContig(unsorted, i));
      }
    }
    Collections.sort(sortMe);
    final NegativeMap translations = new NegativeMap();
    for (OriginalContig contig : sortMe) {
      final long newId = sorted.addContig(contig.mContig);
      for (String attr : unsorted.contigAttributes().keySet()) {
        final String value = unsorted.contigAttribute(contig.mOriginalId, attr);
        if (value != null) {
          sorted.setContigAttribute(newId, attr, value);
        }
      }
      translations.put(contig.mOriginalId, newId);
    }
    final List<TranslatedPath> paths = new ArrayList<>();
    for (int i = 1; i <= unsorted.numberPaths(); ++i) {
      if (!unsorted.pathDeleted(i)) {
        paths.add(new TranslatedPath(unsorted, i, translations));
      }
    }
    Collections.sort(paths);
    for (TranslatedPath path : paths) {
      final long newId = sorted.addPath(path.mPath);
      for (String attr : unsorted.pathAttributes().keySet()) {
        final String value = unsorted.pathAttribute(path.mOriginalId, attr);
        if (value != null) {
          sorted.setPathAttribute(newId, attr, value);
        }
      }
    }
    return sorted;
  }

  /**
   *
   * @param args input graph output directory
   * @throws IOException when you least expect it
   */
  public static void main(String[] args) throws IOException {
    final Graph g = GraphReader.read(new StoreDirProxy(new File(args[0])));
    final File output = new File(args[1]);
    if (output.mkdir()) {
      GraphWriter.write(GraphSorter.sortedGraph(g), new StoreDirProxy(output), "graphsorter", Collections.<UUID>emptySet());
    }
  }
}
