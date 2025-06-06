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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
      final int result = (int) (mOriginalId ^ (mOriginalId >>> 32));
      return 31 * result + (mContig != null ? mContig.hashCode() : 0);
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
      final int result = (int) (mOriginalId ^ (mOriginalId >>> 32));
      return 31 * result + (mPath != null ? mPath.hashCode() : 0);
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
      GraphWriter.write(GraphSorter.sortedGraph(g), new StoreDirProxy(output), "graphsorter", Collections.emptySet());
    }
  }
}
