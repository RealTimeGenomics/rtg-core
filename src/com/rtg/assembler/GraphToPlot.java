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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.LongUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;
import com.rtg.util.store.StoreDirProxy;
/**
 * Read in a graph from our file format and spit out in an alternative visualisation format
 */
public final class GraphToPlot extends LoggedCli {

  static class Node implements Comparable<Node> {
    String mShape;
    final long mContigId;
    String mColor;
    int mLength;
    int mCoverage;
    boolean mTerminal = false;
    Node(long contigId, Graph graph) {
      mContigId = contigId;
      mLength = graph.contigLength(contigId);
      int coverage;
      try  {
        final int kmerFreq = Integer.parseInt(graph.contigAttribute(contigId, GraphKmerAttribute.K_MER_FREQ));
        coverage = kmerFreq / (mLength - graph.contigOverlap());
      } catch (NumberFormatException e) {
        coverage = -1;
      }
      mCoverage = coverage;
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (o == null || getClass() != o.getClass()) {
        return false;
      }
      final Node node = (Node) o;
      return mContigId == node.mContigId;
    }

    @Override
    public int hashCode() {
      return LongUtils.hashCode(mContigId);
    }

    public String output(long id) {
      final StringBuilder sb = new StringBuilder();
      final String coverageString = mCoverage == -1 ? "" : "" + mCoverage;
      final String label = String.format("%s%,d\\n%,d   %s", mContigId > 0 ? "+" : "", mContigId, mLength, coverageString);
      //final String label = (nodeId > 0 ? "+" : "") + nodeId + "\\n" + String.format("%,d", length) + "   " + coverage;
      sb.append("node").append(id).append(" [label=\"").append(label).append("\"");
      if (mColor != null) {
        sb.append(", color=").append(mColor);
      }
      if (mShape != null) {
        sb.append(", shape=").append(mShape);
      }
      if (mTerminal) {
        sb.append(", style=dashed");
      }
      sb.append("];");
      return sb.toString();
    }

    @Override
    public int compareTo(Node o) {
      return Long.compare(mContigId, o.mContigId);
    }
  }

  static class Link implements Comparable<Link> {
    long mFrom;
    long mTo;
    boolean mReverse;
    Link(long from, long to, boolean reverse) {
      mFrom = from;
      mTo = to;
      mReverse = reverse;
    }

    private String output(Graph graph, Map<Long, Long> translated) {
      final StringBuilder sb = new StringBuilder();
      sb.append("node");
      sb.append(translated.get(mFrom));
      sb.append(" -> node");
      sb.append(translated.get(mTo));
      if (graph.pathAttributes().containsKey(GraphKmerAttribute.READ_COUNT)) {
        final int readCount = overlappingReadCount(graph, mFrom, mTo);
        sb.append("[");
        String join  = "";
        if (readCount > 0) {
          sb.append("label=\"").append(readCount).append("\"");
          join = ", ";
        }
        if (mReverse) {
          sb.append(join).append("color=red");
        }
        sb.append("];");

      }
      return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (o == null || getClass() != o.getClass()) {
        return false;
      }

      final Link link = (Link) o;
      return mFrom == link.mFrom && mReverse == link.mReverse && mTo == link.mTo;
    }

    @Override
    public int hashCode() {
      return Utils.pairHashContinuous(LongUtils.hashCode(mFrom), LongUtils.hashCode(mTo), Boolean.valueOf(mReverse).hashCode());
    }

    @Override
    public int compareTo(Link o) {
      return Long.compare(mFrom, o.mFrom) != 0 ? Long.compare(mFrom, o.mFrom) : Long.compare(mTo, o.mTo);
    }
  }

  static final String OUTPUT = "output";
  static final String START = "start";
  static final String WIDTH = "width";

  /**
   * @param args command line
   */
  public static void main(String[] args) {
    new GraphToPlot().mainInit(args, System.out, System.err);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File in = (File) mFlags.getAnonymousValue(0);
    final File output = (File) mFlags.getValue(OUTPUT);
    final Graph graph = GraphReader.read(GraphFactory.KMER, new StoreDirProxy(in));
    final File dotOutput = new File(output, "graph.dot");
    final File pathOutput = new File(output, "long-paths.list");
    final int width = (Integer) mFlags.getValue(WIDTH);
    final long start = (Long) mFlags.getValue(START);

    try (final PrintStream graphStream = new PrintStream(FileUtils.createOutputStream(dotOutput, false));
        final PrintStream pathStream = new PrintStream(FileUtils.createOutputStream(pathOutput, false))) {
      writeDotFile(graph, graphStream, pathStream, width, start);
    }
    return 0;
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);

  }
  protected static void initFlagsLocal(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Produces graphs of the contigs in the specified graph directory");
    flags.registerRequired('o', OUTPUT, File.class, "DIR", "output directory");
    flags.registerRequired(File.class, "DIR", "input graph directory");
    flags.registerRequired('s', START, Long.class, "INT", "produce a .dot file for nodes around this one");
    flags.registerOptional('w', WIDTH, Integer.class, "INT", "maximum distance from the initial node within the .dot", 5);
    flags.setValidator(new Valid());
  }

  @Override
  public String moduleName() {
    return "graph2plot";
  }

  @Override
  public String description() {
    return null;
  }

  private static class Valid implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final File output = (File) flags.getValue(OUTPUT);
      if (!CommonFlags.validateOutputDirectory(output)) {

        return false;
      }
      final File input = (File) flags.getAnonymousValue(0);
      if (!(input.exists() && input.isDirectory())) {
        flags.error("input file should be a directory in the RTG graph file format");
        return false;
      }
      return true;
    }
  }

  static void writeDotFile(Graph graph, PrintStream graphOutput, PrintStream pathStream, int depth, long start) {
    final Set<Node> nodes = collectNodes(graph, depth, start);
    final Set<Long> paths = collectPaths(graph, nodes);
    writeGraph(graph, nodes, paths, graphOutput);
    writeLongPaths(graph, paths, pathStream);
  }

  private static void writeLongPaths(Graph graph, Set<Long> paths, PrintStream out) {
    String joinLine = "";
    for (long path : paths) {
      if (graph.pathLength(path) > 2) {
        out.print(joinLine);
        final String simpleString = graph.pathAttribute(path, GraphKmerAttribute.READ_COUNT);
        final int simpleCount = simpleString == null ? 0 : Integer.parseInt(simpleString);
        final int overlapCount = overlappingReadCount(path, graph);
        out.print("{simpleCount: " + simpleCount  + ", overlapCount: " + overlapCount + ", contigs: [");
        String joinPath = "";
        for (int i = 0; i < graph.pathLength(path); i++)  {
          out.print(joinPath);
          out.print(graph.pathContig(path, i));
          joinPath = ", ";
        }
        out.print("]}");
        joinLine = "," + StringUtils.LS;
      }
    }
    out.println();
  }

  private static void writeGraph(Graph graph, Set<Node> nodes, Set<Long> paths, PrintStream out) {
    out.println("digraph contigGraph {");
    out.println("graph [rankdir=LR, ratio=fill]");
    final Map<Long, Long> translated = new LinkedHashMap<>();
    long i = 0;
    for (Node node : nodes) {
      translated.put(node.mContigId, i);
      out.println(node.output(i));
      i++;
    }

    final Set<Link> links = new TreeSet<>();
    for (long pathId : paths) {
      for (int j = 1; j < graph.pathLength(pathId); j++)  {
        final long first = graph.pathContig(pathId, j - 1);
        final long second = graph.pathContig(pathId, j);
        Link link = null;
        if (translated.containsKey(first) && translated.containsKey(second)) {
          link = new Link(first, second, false);
        } else if (translated.containsKey(first) && translated.containsKey(-second)) {
          link = new Link(first, -second, true);
        } else if (translated.containsKey(-first) && translated.containsKey(second)) {
          link = new Link(-first, second, true);
        }
        if (link != null) {
          links.add(link);
        }
      }
    }
    for (Link link : links) {
      out.println(link.output(graph, translated));
    }
    out.println("}");


  }

  static int overlappingReadCount(long pathId, Graph graph) {
    final long[] path = new long[graph.pathLength(pathId)];
    for (int i = 0; i < path.length; i++) {
      path[i] = graph.pathContig(pathId, i);
    }
    return overlappingReadCount(graph, path);
  }
  static int overlappingReadCount(Graph graph, long... contigs) {
    int count = 0;
    final PathsIterator iterator = graph.paths(contigs[0]);
    long current;
    pathLoop: while ((current = iterator.nextPathId()) != 0) {
      final int index = iterator.contigIndex();
      if (graph.pathLength(current) - index < contigs.length) {
        continue;
      }
      for (int i = 0; i < contigs.length; i++) {
        if (graph.pathContig(current, index + i) != contigs[i]) {
          continue pathLoop;
        }
      }
      final String readCountString = graph.pathAttribute(current, GraphKmerAttribute.READ_COUNT);
      if (readCountString != null) {
        count += Integer.parseInt(readCountString);
      }
    }
    return count;
  }

  private static Set<Long> collectPaths(Graph graph, Set<Node> nodes) {
    final Set<Long> paths = new HashSet<>();
    for (Node node : nodes) {
      final PathsIterator iterator = graph.paths(node.mContigId);
      long pathId;
      while ((pathId = iterator.nextPathId()) != 0) {
        paths.add(pathId);
      }
    }
    return paths;
  }

  static List<Node> nextNodes(long contigId, Graph graph) {
    final List<Node> next = new ArrayList<>();
    final PathsIterator iterator = graph.paths(contigId);
    long pathId;
    while ((pathId = iterator.nextPathId()) != 0) {
      final int index = iterator.contigIndex();
      for (int j = 0; j < graph.pathLength(pathId); j++) {
        if (j == index + 1 || j == index - 1) {
          next.add(new Node(graph.pathContig(pathId, j), graph));
        }
      }

    }
    return next;
  }
  private static Set<Node> collectNodes(Graph graph, int depth, long start) {
    final Set<Node> nodes = new TreeSet<>();
    Set<Node> current = new HashSet<>();
    Set<Node> next = new HashSet<>();
    final Node startNode = new Node(start, graph);
    startNode.mColor = "red";
    startNode.mShape = "box";
    current.add(startNode);
    nodes.add(startNode);
    for (int i = 0; i < depth; i++) {
      for (Node step : current) {
        for (Node linkedContig : nextNodes(step.mContigId, graph)) {
          final Node reverse = new Node(-linkedContig.mContigId, graph);
          if (!nodes.contains(linkedContig) && !nodes.contains(reverse) && !next.contains(reverse)) {
            next.add(linkedContig);
          }
        }
      }
      nodes.addAll(next);
      current = next;
      next = new HashSet<>();
    }
    for (Node node : current) {
      final PathsIterator iterator = graph.paths(node.mContigId);
      long pathId;
      while ((pathId = iterator.nextPathId()) != 0) {
        final int index = iterator.contigIndex();
        for (int j = 0; j < graph.pathLength(pathId); j++) {
          if (j == index + 1 || j == index - 1) {
            final Node linkedContig = new Node(graph.pathContig(pathId, j), graph);
            final Node reverse = new Node(-graph.pathContig(pathId, j), graph);
            if (!nodes.contains(linkedContig) && !nodes.contains(reverse) && !next.contains(reverse)) {
              node.mTerminal = true;
            }
          }
        }
      }
    }
    return nodes;
  }
}
