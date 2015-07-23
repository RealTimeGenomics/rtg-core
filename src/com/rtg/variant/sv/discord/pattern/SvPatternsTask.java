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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;

/**
 */
public class SvPatternsTask extends ParamsTask<BreakpointPatternParams, NoStatistics> implements Closeable {
  //private static Pattern PATTERN = Pattern.compile("\\s+");
  private final int mSameDistance;
//  private final int mMaxScan;
  private final int mMinDepth;

  private OutputStream mOutput;
  private final List<DiscordBedRecord> mBedRecords = new ArrayList<>();
  private final DeletionOverlapFilter mDeletions = new DeletionOverlapFilter();
  private static final String FILE_NAME = "sv_patterns.bed";
  private static final String SV_OUTPUT_VERSION = "0.1";

  boolean filterRecord(VcfRecord rec) {
    if (rec.isFiltered()) {
      return true;
    }
    final Collection<String> depth = rec.getInfo().get("DP");
    for (String d : depth) {
      try {
        if (Integer.parseInt(d) < mMinDepth) {
          return true;
        }
      } catch (NumberFormatException e) {
        return true;
      }
    }
    return false;
  }
  BreakpointStore loadBreakpoints(RegionRestriction region,  List<File> files) throws IOException {
    final BreakpointStore store = new BreakpointStore();
    for (File f : files) {

      final VcfReader vcf = VcfReader.openVcfReader(f, region);
      while (vcf.hasNext()) {
        final VcfRecord record =  vcf.next();
        if (!filterRecord(record)) {
          store.add(new VcfBreakpoint(record));
        }
      }
    }
    return store;
  }

  void analyseBreakpoint(VcfBreakpoint breakpoint, BreakpointStore store) {
    final String local = breakpoint.getLocalChr();
    final String remote = breakpoint.getRemoteChr();
    if (local.equals(remote) && breakpoint.isLocalUp() && !breakpoint.isRemoteUp()) {
      potentialDelete(breakpoint);
    }
    if (local.equals(remote) && breakpoint.isLocalUp() && breakpoint.isRemoteUp() && breakpoint.getRemotePos() > breakpoint.getLocalPos()) {
      potentialInversion(breakpoint, store);
    }
    if (breakpoint.isLocalUp()) {
      potentialCopy(breakpoint, store);
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
  private void potentialDelete(VcfBreakpoint breakpoint) {
    if (breakpoint.getLocalPos() < breakpoint.getRemotePos()) {
      final DiscordBedRecord bed = new DiscordBedRecord(breakpoint.getLocalChr(), breakpoint.getLocalPos(), breakpoint.getRemotePos(),
          "deletion",
          "" + breakpoint.getDepth());
      write(bed);
    }
  }

  void write(DiscordBedRecord bed) {
    if (bed.getAnnotations()[0].equals("deletion")) {
      mDeletions.add(bed);
    } else {
      mBedRecords.add(bed);
    }
  }

  SvPatternsTask(BreakpointPatternParams params, OutputStream output) {
    super(params, output, new NoStatistics(), null);
    mSameDistance = params.sameDistance();
    //mMaxScan = params.fragmentLength();
    mMinDepth = params.minDepth();
  }

  @Override
  protected void exec() throws IOException {
    final File bedFile = new File(mParams.directory(), FILE_NAME + FileUtils.GZ_SUFFIX);
    setOutput(FileUtils.createOutputStream(bedFile, true));
    try {
      writeHeader(mOutput);
      final RegionRestriction region = mParams.region();
      final List<File> files = mParams.files();
      final BreakpointStore store = loadBreakpoints(region, files);
      processStore(store);
    } finally {
      close();
    }

    final File index = new File(bedFile.getPath() + TabixIndexer.TABIX_EXTENSION);
    final TabixIndexer tabixIndexer = new TabixIndexer(bedFile, index);
    try {
      tabixIndexer.saveBedIndex();
    } catch (UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + bedFile + ": " + e.getMessage());
    }
  }

  void processStore(BreakpointStore store) throws IOException {
    for (VcfBreakpoint breakpoint : store) {
      analyseBreakpoint(breakpoint, store);
    }
  }

  private void writeHeader(OutputStream output) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("#Version ").append(Environment.getVersion()).append(", SV Patterns output ").append(SvPatternsTask.SV_OUTPUT_VERSION).append(LS);
    if (CommandLine.getCommandLine() != null) {
      sb.append("#CL" + TAB).append(CommandLine.getCommandLine()).append(LS);
    }
    sb.append("#RUN-ID" + TAB).append(CommandLine.getRunId()).append(LS);
    sb.append("#");
    sb.append("sequence" + TAB + "start" + TAB + "end" + TAB);
    sb.append("description" + TAB + "count");
    sb.append(LS);
    output.write(sb.toString().getBytes());
  }
  void setOutput(OutputStream out) {
    mOutput = out;
  }
  @Override
  public void close() throws IOException {
    mBedRecords.addAll(mDeletions.nonOverlapping());
    Collections.sort(mBedRecords, new SmartBedWriter.BedPositionalComparator());
    for (DiscordBedRecord b : mBedRecords) {
      mOutput.write(b.toString().getBytes());
      ByteUtils.writeLn(mOutput);
    }
    mOutput.close();
  }

}
