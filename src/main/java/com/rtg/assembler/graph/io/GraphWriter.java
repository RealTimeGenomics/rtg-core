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

package com.rtg.assembler.graph.io;

import static com.rtg.util.StringUtils.TAB;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.UUID;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Graph;
import com.rtg.util.MD5Utils;
import com.rtg.util.store.StoreDirectory;
import com.rtg.mode.DNARange;

/**
 * Write a graph to a directory in version 0.0 of the standard file format.
 */
@TestClass({"com.rtg.assembler.graph.io.GraphWriterTest", "com.rtg.assembler.graph.io.CommonIOTest"})
public final class GraphWriter {

  private final boolean mIncludeDeletions;

  /**
   * Write a standard version 0.0 graph into the directory.
   * @param graph to be written.
   * @param dir the directory where everything is to be written.
   * @param command command used to generate the graph.
   * @param inputUids <code>UUID</code>s of inputs to the program that created the graph.
   * @throws IOException whenever.
   */
  public static void write(final Graph graph, final StoreDirectory dir, final String command, final Set<UUID> inputUids) throws IOException {
    new GraphWriter(graph, dir, Integer.MAX_VALUE, false).write(command, inputUids);
  }

  /**
   * Write a version 0.0 graph preserving information about deleted contigs/paths
   * @param graph to be written.
   * @param dir the directory where everything is to be written.
   * @param command command used to generate the graph.
   * @param inputUids <code>UUID</code>s of inputs to the program that created the graph.
   * @throws IOException whenever.
   */
  public static void writeWithDeleted(final Graph graph, final StoreDirectory dir, final String command, final Set<UUID> inputUids) throws IOException {
    new GraphWriter(graph, dir, Integer.MAX_VALUE, true).write(command, inputUids);
  }

  private final Graph mGraph;

  private final StoreDirectory mDir;

  private final int mMaxLength;

  //Only use internally and for testing
  GraphWriter(final Graph graph, final StoreDirectory dir, final int maxLength, boolean deleted) {
    mGraph = graph;
    mDir = dir;
    mMaxLength = maxLength;
    mIncludeDeletions = deleted;
  }

  //Exposed only for testing
  void write(final String command, final Set<UUID> inputUids) throws IOException {
    writeContig();
    writePath();
    writeHeader(mGraph, command, new Date(), UUID.randomUUID(), inputUids);
    writeDigest();
  }

  private void writeContig() {
    try (final PrintStream ps = new PrintStream(new BoundedStreams(mDir, mMaxLength, "contig.", ".fa"))) {
      final Set<String> contigAttributes = new TreeSet<>(mGraph.contigAttributes().keySet());
      for (long l = 1; l <= mGraph.numberContigs(); ++l) {
        if (!mIncludeDeletions && mGraph.contigDeleted(l)) {
          continue;
        }
        ps.print(">" + l);
        for (final String key : contigAttributes) {
          if (mIncludeDeletions || !"deleted".equals(key)) {
            final String attr = mGraph.contigAttribute(l, key);
            if (attr != null) {
              ps.print(TAB + key + "=" + attr);
            }
          }
        }
        ps.println();
        final int length = mGraph.contigLength(l);
        for (int i = 0; i < length; ++i) {
          if (i > 0 && (i % 80) == 0) {
            ps.println();
          }
          final byte nt = mGraph.nt(l, i);
          ps.print(DNARange.RANGE.toChar(nt));
        }
        ps.println();
        ps.flush();
      }
    }
  }
  private void writePath() {
    try (final PrintStream ps = new PrintStream(new BoundedStreams(mDir, mMaxLength, "path.", ".tsv"))) {
      final Set<String> pathAttributes = new TreeSet<>(mGraph.pathAttributes().keySet());
      for (long l = 1; l <= mGraph.numberPaths(); ++l) {
        if (!mIncludeDeletions && mGraph.pathDeleted(l)) {
          continue;
        }
        ps.print("path" + TAB + l);
        for (final String key : pathAttributes) {
          if (mIncludeDeletions || !"deleted".equals(key)) {
            final String attr = mGraph.pathAttribute(l, key);
            if (attr != null) {
              ps.print(TAB + key + "=" + attr);
            }
          }
        }
        ps.println();
        ps.print("contigs");
        final int length = mGraph.pathLength(l);
        for (int i = 0; i < length; ++i) {
          final long contig = mGraph.pathContig(l, i);
          ps.print(TAB + contig);
        }
        ps.println();
        ps.flush();
      }
    }
  }

  private static void date(final Date date, final PrintStream ps) {
    ps.print(dateFormat().format(date));
  }

  static SimpleDateFormat dateFormat() {
    return new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
  }

  void writeHeader(final Graph graph, final String command, final Date date, final UUID uuid, final Set<UUID> inputGuids) throws IOException {
    try (final PrintStream ps = new PrintStream(mDir.child("header.tsv").outputStream())) {
      ps.println("version" + TAB + "0.0");

      ps.print("date" + TAB);
      date(date, ps);
      ps.println();

      ps.println("command" + TAB + command);

      ps.println("guid" + TAB + uuid.toString());
      ps.println("contigOverlap" + TAB + graph.contigOverlap());

      for (final UUID inputGuid : new TreeSet<>(inputGuids)) {
        ps.println("inputguid" + TAB + inputGuid);
      }

      writeAttributeKeys(ps, "contig", graph.contigAttributes());
      writeAttributeKeys(ps, "path", graph.pathAttributes());
    }
  }

  static void writeAttributeKeys(final PrintStream ps, final String key, final Map<String, String> attributes) {
    for (final String entry : new TreeSet<>(attributes.keySet())) {
      ps.println(key + TAB + entry + TAB + attributes.get(entry));
    }
  }

  private void writeDigest() throws IOException {
    final SortedSet<String> children = mDir.children(); //This line must go before the next to prevent the manifest itself being included.
    try (final PrintStream ps = new PrintStream(mDir.child("manifest.txt").outputStream())) {
      for (final String name : children) {
        try (final InputStream is = mDir.child(name).inputStream()) {
          final String md5 = MD5Utils.md5(is);
          ps.println(md5 + " " + name);
        }
      }
    }
  }
}
