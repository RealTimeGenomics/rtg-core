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

import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.io.LogStream;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public class FilterPaths extends LoggedCli {

  static final String READ_COUNT = "read-count";
  static final String MIN_PATH = "min-path";
  static final String INCLUDE_DELETED = "include-deleted";


  /**
   * @param args command line
   */
  public static void main(String[] args) {
    final FilterPaths filterPaths = new FilterPaths();
    CommandLine.setCommandArgs(args);
    filterPaths.mainInit(args, System.out, System.err);
  }

  @Override
  public String moduleName() {
    return "filterpaths";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);

  }
  static void initCommonFlags(CFlags flags) {
    flags.registerOptional('r', READ_COUNT, Integer.class, "Int", "minimum reads required to overlap a link").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('p', MIN_PATH, Integer.class, "Int", "minimum reads required in long paths").setCategory(SENSITIVITY_TUNING);
  }
  protected static void initFlagsLocal(CFlags flags) {
    flags.setDescription("Attempts to reduce the complexity of a graph by removing invalid paths");
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    initCommonFlags(flags);
    flags.registerOptional('d', INCLUDE_DELETED, "set if you should output nodes that are deleted").setCategory(INPUT_OUTPUT);
    flags.registerRequired(File.class, CommonFlags.DIR, "input graph directory").setCategory(INPUT_OUTPUT);
    flags.setValidator(new Valid());
  }

  private static class Valid implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      final File input = (File) flags.getAnonymousValue(0);
      if (!(input.exists() && input.isDirectory())) {
        flags.setParseMessage("Input file should be a directory in the RTG graph file format");
        return false;
      }
      return true;
    }
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File in = (File) mFlags.getAnonymousValue(0);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final StoreDirProxy graphDir = new StoreDirProxy(in);
    final Set<UUID> sourceIds = new HashSet<>();
    sourceIds.add(GraphReader.getUUID(graphDir));
    final Graph graph = GraphReader.read(GraphFactory.KMER, graphDir);
    final MutableGraph mutable = (MutableGraph) graph;

    if (mFlags.isSet(READ_COUNT)) {
      final int readCount = (Integer) mFlags.getValue(READ_COUNT);
      improveSingle(mutable, readCount);
    }
    if (mFlags.isSet(MIN_PATH)) {
      final int minPath = (Integer) mFlags.getValue(MIN_PATH);
      improveMultiple(mutable, minPath);
    }

    if (mFlags.isSet(INCLUDE_DELETED)) {
      GraphWriter.writeWithDeleted(mutable, new StoreDirProxy(output), CommandLine.getCommandLine(), sourceIds);
    } else {
      GraphWriter.write(mutable, new StoreDirProxy(output), CommandLine.getCommandLine(), sourceIds);

    }
    return 0;
  }

  static void improveSingle(MutableGraph mutable, int threshold) {
    for (int i = 1; i <= mutable.numberPaths(); ++i) {
      if (GraphToPlot.overlappingReadCount(i, mutable) < threshold) {
        mutable.deletePath(i);
      }
    }
  }
  static int readCount(long pathId, Graph graph) {
    final String attr = graph.pathAttribute(pathId, GraphKmerAttribute.READ_COUNT);
    if (attr == null) {
      return 0;
    } else {
      return Integer.parseInt(attr);
    }
  }
  static void improveMultiple(MutableGraph graph, int threshold) {
    for (long i = 1; i < graph.numberPaths(); ++i) {
      if (graph.pathDeleted(i) || graph.pathLength(i) != 2) {
        continue;
      }
      final List<Long> unique = biDirectional(graph.pathContig(i, 0), graph.pathContig(i, 1), graph, threshold);
      if (unique.size() > 2) {
        final long otherEnd = GraphMap.findPath(unique.subList(unique.size() - 2, unique.size()), graph, true);
        if (otherEnd != 0) {
          graph.deletePath(otherEnd);
        }
        graph.deletePath(i);
      }
    }

  }
  static List<Long> biDirectional(long first, long second, Graph graph, int threshold) {
    final List<Long> forward = unambiguousPath(first, second, graph, threshold);
    for (int end = forward.size() - 1; end >= 1; --end) {
      final List<Long> reverse = unambiguousPath(-forward.get(end), -forward.get(end - 1), graph, threshold);
      if (reverse.size() > end && reverse.get(end) == -first) {
        return forward.subList(0, end + 1);
      }
    }
    return Collections.emptyList();
  }
  static List<Long> unambiguousPath(long first, long second, Graph graph, int threshold) {
    int failPosition = Integer.MAX_VALUE;
    final List<Integer> counts = new ArrayList<>();
    final List<Long> path = new ArrayList<>();
    path.add(first);
    path.add(second);
    counts.add(0);
    counts.add(0);
    final PathsIterator iterator = graph.paths(first);
    long current;
    while ((current = iterator.nextPathId()) != 0) {
      final int index = iterator.contigIndex();
      if (graph.pathLength(current) > index + 1 && graph.pathContig(current, index + 1) != second) {
        continue;
      }
      for (int contig = index, pathPos = 0; contig < graph.pathLength(current); ++contig, ++pathPos) {
        final int readCount = readCount(current, graph);
        if (pathPos > path.size() - 1) {
          path.add(graph.pathContig(current, contig));
          counts.add(readCount);
        } else if (path.get(pathPos) != graph.pathContig(current, contig)) {
          failPosition = Math.min(failPosition, pathPos);
          break;
        } else {
          final int newCount = counts.get(pathPos) + readCount;
          counts.set(pathPos, newCount);
        }
      }
    }
    for (int i = 0; i < failPosition && i < path.size(); ++i) {
      if (counts.get(i) < threshold) {
        return path.subList(0, i);
      }
    }
    return path.subList(0, Math.min(failPosition, path.size()));
  }

}
