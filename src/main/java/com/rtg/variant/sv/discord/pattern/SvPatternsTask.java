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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.NavigableSet;
import java.util.SortedSet;

import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.ByteUtils;
import com.rtg.util.Environment;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.sv.discord.DiscordBedRecord;
import com.rtg.variant.sv.discord.SmartBedWriter;
import com.rtg.vcf.AltVariantTypeFilter;
import com.rtg.vcf.AssertVcfSorted;
import com.rtg.vcf.PassOnlyFilter;
import com.rtg.vcf.VariantType;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfFilterIterator;
import com.rtg.vcf.VcfInfoFilter;
import com.rtg.vcf.VcfIterator;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfUtils;

/**
 */
public class SvPatternsTask extends ParamsTask<BreakpointPatternParams, NoStatistics> implements Closeable {

  private static final String FILE_NAME = "sv_patterns.bed";
  private static final String SV_OUTPUT_VERSION = "0.1";

  private final int mSameDistance;
  private final int mMinDepth;

  private OutputStream mOutput;
  private final List<DiscordBedRecord> mBedRecords = new ArrayList<>();
  private final DeletionOverlapFilter mDeletions = new DeletionOverlapFilter();

  SvPatternsTask(BreakpointPatternParams params, OutputStream output) {
    super(params, output, new NoStatistics(), null);
    mSameDistance = params.sameDistance();
    mMinDepth = params.minDepth();
  }

  @Override
  protected void exec() throws IOException {
    final File bedFile = new File(mParams.directory(), FILE_NAME + FileUtils.GZ_SUFFIX);
    setOutput(FileUtils.createOutputStream(bedFile));
    try {
      writeHeader(mOutput);
      final BreakpointStore store = loadBreakpoints(mParams.region(), mParams.files());
      for (VcfBreakpoint breakpoint : store) {
        analyseBreakpoint(breakpoint, store);
      }
    } finally {
      close();
    }

    try {
      new TabixIndexer(bedFile).saveBedIndex();
    } catch (UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + bedFile + ": " + e.getMessage());
    }
  }

  void setOutput(OutputStream out) {
    mOutput = out;
  }

  BreakpointStore loadBreakpoints(RegionRestriction region, List<File> files) throws IOException {
    final List<VcfFilter> filters = new ArrayList<>();
    filters.add(new AssertVcfSorted());
    filters.add(new PassOnlyFilter());
    filters.add(new AltVariantTypeFilter(EnumSet.of(VariantType.SV_BREAKEND)));
    filters.add(new VcfInfoFilter.MinMaxIntFilter(null, null, mMinDepth, Integer.MAX_VALUE, VcfUtils.INFO_COMBINED_DEPTH));

    final BreakpointStore store = new BreakpointStore();
    for (File f : files) {
      try (final VcfIterator vcf = new VcfFilterIterator(VcfReader.openVcfReader(f, region), filters)) {
        vcf.forEach(rec -> store.add(new VcfBreakpoint(rec)));
      }
    }
    return store;
  }

  void analyseBreakpoint(VcfBreakpoint breakpoint, BreakpointStore store) {
    final String local = breakpoint.getLocalChr();
    final String remote = breakpoint.getRemoteChr();
    if (local.equals(remote)) {
      if (breakpoint.isLocalUp() && !breakpoint.isRemoteUp() && breakpoint.getLocalPos() < breakpoint.getRemotePos()) {
        potentialDeleteUp(breakpoint);
      } else if (!breakpoint.isLocalUp() && breakpoint.isRemoteUp() && breakpoint.getLocalPos() > breakpoint.getRemotePos()) {
        potentialDeleteDown(breakpoint);
      }
      if (breakpoint.isLocalUp() && breakpoint.isRemoteUp() && breakpoint.getRemotePos() > breakpoint.getLocalPos()) {
        potentialInversion(breakpoint, store);
      }
    }
    if (breakpoint.isLocalUp()) {
      potentialCopy(breakpoint, store);
    }
  }

  // Assume reference span of 1, which needs to be accounted for in start positions (effectively to remove padding base)
  // TODO, check other event types for potential adjustment.
  private static final int REF_SPAN_ADJUST = 1;

  private void potentialDeleteUp(VcfBreakpoint breakpoint) {
    final DiscordBedRecord bed = new DiscordBedRecord(breakpoint.getLocalChr(), breakpoint.getLocalPos() + REF_SPAN_ADJUST, breakpoint.getRemotePos(),
      "deletion", "" + breakpoint.getDepth());
    mDeletions.add(bed);
  }
  private void potentialDeleteDown(VcfBreakpoint breakpoint) {
    final DiscordBedRecord bed = new DiscordBedRecord(breakpoint.getLocalChr(), breakpoint.getRemotePos() + REF_SPAN_ADJUST, breakpoint.getLocalPos(),
      "deletion", "" + breakpoint.getDepth());
    mDeletions.add(bed);
  }

  private void potentialInversion(VcfBreakpoint breakpoint, BreakpointStore store) {
    assert breakpoint.isLocalUp() : "We're only handling the first edge of breakpoints to prevent duplication";
    final BreakpointStore.PositionMap similar = store.getMap().get(breakpoint.getLocalChr()).get(breakpoint.getRemoteChr());
    final VcfBreakpoint start = firstSame(breakpoint, similar);
    // Larger inversions
    final NavigableSet<VcfBreakpoint> tail = similar.tailSet(start, true);
    for (VcfBreakpoint b : tail) {
      if (b == breakpoint) {
        continue;
      }
      if (b.getLocalPos() > breakpoint.getLocalPos() + mSameDistance) {
        break;
      }
      if (!b.isLocalUp() && !b.isRemoteUp() && (Math.abs(b.getRemotePos() - breakpoint.getRemotePos()) < mSameDistance)) {
        final DiscordBedRecord bed = new DiscordBedRecord(breakpoint.getLocalChr(), breakpoint.getLocalPos(), breakpoint.getRemotePos(),
            "inversion",
            "" + (breakpoint.getDepth() + b.getDepth()));
        write(bed);
      }
    }
  }

  private void potentialCopy(VcfBreakpoint breakpoint, BreakpointStore store) {
    assert breakpoint.isLocalUp() : "We're only handling the first edge of breakpoints to prevent duplication";
    // map of breakpoints across the same chromosomes
    final BreakpointStore.PositionMap similar = store.getMap().get(breakpoint.getLocalChr()).get(breakpoint.getRemoteChr());
    final SortedSet<VcfBreakpoint> tail = similar.tailSet(firstSame(breakpoint, similar), true);
    for (VcfBreakpoint other: tail) {
      if (other.getLocalPos() > breakpoint.getLocalPos() + mSameDistance) {
        break;
      }
      // other is pointing in the other direction at both sides of the breakpoint
      if (!other.isLocalUp() && other.isRemoteUp() != breakpoint.isRemoteUp()) {
        //
        if (breakpoint.isRemoteUp() && other.getRemotePos() < breakpoint.getRemotePos()
            || !breakpoint.isRemoteUp() && other.getRemotePos() > breakpoint.getRemotePos()) {
          // We might be legit
          final int start = Math.min(breakpoint.getLocalPos(), other.getLocalPos());
          final int end = Math.max(breakpoint.getLocalPos(), other.getLocalPos());
          final DiscordBedRecord bed = new DiscordBedRecord(breakpoint.getLocalChr(), start, end,
              "inserted_copy" + ":" + breakpoint.getRemoteChr() + ":" + breakpoint.getRemotePos() + "-" + other.getRemotePos(),
              "" + (breakpoint.getDepth() + other.getDepth()));
          write(bed);

        }
      }
    }
  }

  private VcfBreakpoint firstSame(VcfBreakpoint breakpoint, BreakpointStore.PositionMap similar) {
    final NavigableSet<VcfBreakpoint> head = similar.headSet(breakpoint, false);
    VcfBreakpoint start = breakpoint;
    for (VcfBreakpoint b : head.descendingSet()) {
      if (b.getLocalPos() < breakpoint.getLocalPos() - mSameDistance) {
        break;
      }
      start = b;
    }
    return start;
  }

  void write(DiscordBedRecord bed) {
    mBedRecords.add(bed);
  }

  private void writeHeader(OutputStream output) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("#Version ").append(Environment.getVersion()).append(", SV Patterns output ").append(SvPatternsTask.SV_OUTPUT_VERSION).append(LS);
    if (CommandLine.getCommandLine() != null) {
      sb.append("#CL" + TAB).append(CommandLine.getCommandLine()).append(LS);
    }
    sb.append("#");
    sb.append("sequence" + TAB + "start" + TAB + "end" + TAB);
    sb.append("description" + TAB + "count");
    sb.append(LS);
    output.write(sb.toString().getBytes());
  }

  @Override
  public void close() throws IOException {
    mBedRecords.addAll(mDeletions.nonOverlapping());
    mBedRecords.sort(new SmartBedWriter.BedPositionalComparator());
    for (DiscordBedRecord b : mBedRecords) {
      mOutput.write(b.toString().getBytes());
      ByteUtils.writeLn(mOutput);
    }
    mOutput.close();
  }

}
