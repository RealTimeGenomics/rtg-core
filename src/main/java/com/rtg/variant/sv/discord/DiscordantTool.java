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
package com.rtg.variant.sv.discord;


import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.sam.BamIndexer;
import com.rtg.sam.SamIteratorTask;
import com.rtg.sam.SamUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.CompareHelper;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.sv.ReadGroupStats;
import com.rtg.vcf.DefaultVcfWriter;
import com.rtg.vcf.ReorderingVcfWriter;
import com.rtg.vcf.VcfRecord;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Primary task that process SAM records to identify putative structural variant breakends.
 */
public class DiscordantTool extends SamIteratorTask<DiscordantToolParams, DiscordantToolStatistics> {

  /** Specifies what type of BAM output to produce */
  public enum BamType {
    /** Do not write BAM */
    NONE,
    /** Write one BAM file per discordant read set */
    PER_CLUSTER,
    /** Write one BAM file containing all discordant read sets */
    COMBINED
  }

  static final String OUTPUT_FILENAME = "discordant_pairs.vcf";
  static final String DEBUG_FILENAME = "debug_discordant_pairs.tsv";
  static final String BED_FILENAME = "discordant_pairs.bed";
  static final String BAM_FILENAME_PREFIX = "discordant";
  static final String DT_OUTPUT_VERSION = "1";

  // If set, debug mode will also output a separate debug file for every new discordant read.
  // This allows plotting breakpoint constraint evolution.
  // Warning, this buffers all read sets, so only use for smallish regions.
  private static final boolean DEBUG_PER_RECORD = GlobalFlags.isSet(CoreGlobalFlags.SV_DISCORD_DEBUG_PER_RECORD);

  private final List<DiscordantReadSet> mWritten = DEBUG_PER_RECORD ? new ArrayList<>() : null;

  private int mMaxGap;
  protected final SortedSet<DiscordantReadSet> mReadSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());

  final Map<String, MachineOrientation> mMachineOrientations = new HashMap<>();

  private final VcfDiscordantOutputFormatter mFormatter;
  private final DebugDiscordantOutputFormatter mDebugFormatter;
  private final BedDiscordantOutputFormatter mBedFormatter;
  private ReorderingDebugOutput mDebugReordering;
  private final TreeSet<DiscordantReadSet> mOutputBuffer = new TreeSet<>(new BreakpointPositionComparator());
  private final BamType mBamType;

  private OutputStream mOut = null;
  private OutputStream mDebugOutput = null;
  private OutputStream mBedOutput = null;
  private SmartBedWriter mBedWriter = null;
  private ReorderingVcfWriter mVcfWriter = null;

  private SAMFileHeader mSamHeader = null;
  private long mTotalDiscordantRecords = 0;
  private SAMFileWriter mBamWriter = null;
  private File mBamFile = null;

  protected DiscordantTool(DiscordantToolParams params, OutputStream defaultOutput) throws IOException {
    super(params, defaultOutput, new DiscordantToolStatistics(params.directory()), params.filterParams());
    mBamType = params.bamOutput();
    mFormatter = new VcfDiscordantOutputFormatter(mGenomeSequences);
    mBedFormatter = new BedDiscordantOutputFormatter();
    mDebugFormatter =  new DebugDiscordantOutputFormatter();
  }

  private BreakpointConstraint getConstraint(DiscordantReadSet rs) {
    return rs.getIntersection() == null ? rs.getUnion() : rs.getIntersection();
  }
  private int readSetPosition(DiscordantReadSet rs) {
    final BreakpointConstraint bc = getConstraint(rs);
    return Math.min(mTemplateLength, Math.max(0, bc.position().position()));
  }

  @Override
  protected void finalPostFlush() throws IOException {
    flush(0, -1);
    if (mBedWriter != null) {
      mBedWriter.close();
    }
    mVcfWriter.close();
    closeAndIndexBam();
  }

  private void closeAndIndexBam() throws IOException {
    if (mBamWriter != null) {
      mBamWriter.close();
      mBamWriter = null;
      try {
        BamIndexer.saveBamIndex(mBamFile, BamIndexer.indexFileName(mBamFile));
      } catch (final UnindexableDataException e) {
        Diagnostic.warning("Cannot produce index for: " + mBamFile + ": " + e.getMessage());
      }
    }
  }

  @Override
  public int flush(int start, int last) throws IOException {
    // we don't really care about the start... flush anything more than a fragment length before last.
    int res = start;
    if (last == -1) {
      for (final DiscordantReadSet drs : mReadSets) {
        flush(drs);
        res = drs.flushPosition();
      }
      mReadSets.clear();
    } else {
      final Iterator<DiscordantReadSet> it = mReadSets.iterator();
      while (it.hasNext()) {
        final DiscordantReadSet drs = it.next();
        if (drs.flushPosition() >= last) {
          break;
        }
        it.remove();
        flush(drs);
        res = drs.flushPosition();
      }
    }
    output();
    return res;
  }

  private void flush(DiscordantReadSet drs) {
    if (drs.getCounts() >= mParams.minBreakpointDepth() && (!mParams.intersectionOnly() || drs.getIntersection() != null)) {
      mOutputBuffer.add(drs);
    }
  }
  private void output() throws IOException {
    while (!mOutputBuffer.isEmpty()) {
      final DiscordantReadSet first = mOutputBuffer.pollFirst();
      writeReadSet(first, -1, -1);
    }
  }

  private void writeReadSet(DiscordantReadSet drs, int coverage, double ambiguous) throws IOException {
    if (mParams.debugOutput()) {
      mDebugReordering.addRecord(drs);
    }
    if (DEBUG_PER_RECORD) {
      mWritten.add(drs); // Buffers forever, only use on small datasets
    }
    if (mBamType != BamType.NONE) {
      writeReadSetBam(drs);
    }
    final VcfRecord rec = mFormatter.vcfRecord(drs, coverage, ambiguous);
    mVcfWriter.addRecord(rec);
    mStatistics.tallyDiscordantReadSet(drs);
    if (mBedOutput != null) {
      mBedWriter.addRecord(mBedFormatter.format(drs));
    }
  }

  private void writeReadSetBam(DiscordantReadSet drs) throws IOException {
    if (mBamType == BamType.PER_CLUSTER && mBamWriter != null) {
      closeAndIndexBam();
    }
    if (mBamWriter == null) {
      final SAMFileHeader header = mSamHeader.clone();
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      if (mBamType == BamType.PER_CLUSTER) {
        header.addComment(mDebugFormatter.format(drs).trim());
        mBamFile = mParams.file(BAM_FILENAME_PREFIX + "." + drs.getSequenceName() + ":" + drs.unionPosition() + SamUtils.BAM_SUFFIX);
      } else {
        mBamFile = mParams.file(BAM_FILENAME_PREFIX + SamUtils.BAM_SUFFIX);
      }
      mBamWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, FileUtils.createOutputStream(mBamFile), false);
    }
    for (final SAMRecord rec : drs.getRecords()) {
      mBamWriter.addAlignment(rec);
    }
    drs.getRecords().clear();
  }

  @Override
  public void init(SAMFileHeader header) throws IOException {
    mSamHeader = header;
    int maxMaxGap = 0;

    final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
    if (readGroups == null || readGroups.isEmpty()) {
      throw new NoTalkbackSlimException("SAM input does not contain read groups");
    }
    final Set<String> sampleNames = new HashSet<>(readGroups.size());
    for (final SAMReadGroupRecord srgr : readGroups) {
      if (srgr.getPlatform() == null) {
        throw new NoTalkbackSlimException("Read group does not contain a platform field: " + srgr.getSAMString());
      }
      if (srgr.getSample() == null) {
        throw new NoTalkbackSlimException("Read group does not contain a sample field: " + srgr.getSAMString());
      }
      sampleNames.add(srgr.getSample());
      final String platform = srgr.getPlatform().toUpperCase(Locale.ROOT);
      if (MachineType.ILLUMINA_PE.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.ILLUMINA_PE.orientation());
      } else if (MachineType.COMPLETE_GENOMICS.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.COMPLETE_GENOMICS.orientation());
      } else if (MachineType.COMPLETE_GENOMICS_2.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.COMPLETE_GENOMICS_2.orientation());
      } else {
        throw new NoTalkbackSlimException("Unsupported platform: " + srgr.getPlatform());
      }
      final ReadGroupStats rgs = mParams.readGroupStatistics().get(srgr.getId());
      if (rgs == null) {
        throw new NoTalkbackSlimException("Read group " + srgr.getId() + " did not have corresponding rgstats data supplied");
      }

      rgs.setNumDeviations(mParams.numDeviations());
      final double dev = rgs.concordantDeviation();
      Diagnostic.userLog("Read group " + srgr.getId() + " concordant deviation +/-" + MathUtils.round(dev)
        + ", fragment (" + (int) (rgs.fragmentMean() - dev) + "," + (int) rgs.fragmentMean() + "," + (int) (rgs.fragmentMean() + dev) + ")"
        + ", gap (" + (int) (rgs.gapMean() - dev) + "," + (int) rgs.gapMean() + "," + (int) (rgs.gapMean() + dev) + ")");
      maxMaxGap = Math.max(maxMaxGap, rgs.gapMax());
    }
    if (!mParams.multisample() && sampleNames.size() > 1) {
      throw new NoTalkbackSlimException("Input read groups contain multiple samples: " + sampleNames);
    }
    final String vcfSample = sampleNames.size() == 1 ? sampleNames.iterator().next() : "SAMPLE";
    mMaxGap = maxMaxGap;
    if (mDebugOutput != null) {
      mDebugOutput.write(mDebugFormatter.header().getBytes());
      mDebugReordering = new ReorderingDebugOutput(mDebugFormatter, mDebugOutput, 5 * mMaxGap);
    }
    if (mBedOutput != null) {
      mBedOutput.write(mBedFormatter.header().getBytes()); // Maybe move this into the bedwriter
      mBedWriter = new SmartBedWriter(mBedOutput);
    }

    mVcfWriter = new ReorderingVcfWriter(new DefaultVcfWriter(mFormatter.header(header, vcfSample, false, false), mOut));
  }

  @Override
  protected boolean processRecord(final SAMRecord rec) throws IOException {
    if (!rec.getReferenceName().equals(mTemplateName)) {
      mTemplateName = rec.getReferenceName();
      flush(mPreviousStart, -1);
      mTemplateLength = getSequenceLength(mGenomeSequences, mTemplateNameMap, mTemplateName);
    }
    if (rec.getMateUnmappedFlag() || "*".equals(rec.getMateReferenceName())) {
      return true;
    }
    if (SamUtils.uniquelyMapped(rec)) {
      final SAMReadGroupRecord rg = rec.getReadGroup();
      if (rg == null) {
        throw new NoTalkbackSlimException("Input SAM record did not contain a read group: " + rec.getSAMString());
      }
      final String id = rg.getId();
      final MachineOrientation mo = mMachineOrientations.get(id);
      final ReadGroupStats rgs = mParams.readGroupStatistics().get(id);

      // BWA sometimes likes to align tiny fractions of a read (e.g. a 20-mer homopolymer), which are often bollocks.
      // Let's require at least a third of the read bases to be aligned.
      if ((rec.getAlignmentEnd() - rec.getAlignmentStart()) >= rgs.meanLength() / 3) {

        final BreakpointConstraint constraint = new BreakpointConstraint(rec, mo, rgs, mParams.overlapFraction());
        if (!constraint.isConcordant(rgs)) {
          processConstraint(constraint, mReadSets, mTemplateName, mMaxGap, mBamType == BamType.NONE ? null : rec);
          ++mTotalDiscordantRecords;
          if (DEBUG_PER_RECORD) {
            dumpAllReadSets(constraint);
          }
        }
      }
    }
    flush(0, rec.getAlignmentStart());
    return true;
  }

  private void dumpAllReadSets(BreakpointConstraint constraint) throws IOException {
    final File dumpFile = mParams.outFile("debug_discordant_pairs_" + mTotalDiscordantRecords + ".tsv");
    try (OutputStream debugOutput = FileUtils.createOutputStream(dumpFile)) {
      debugOutput.write(mDebugFormatter.header().getBytes());
      try (final ReorderingDebugOutput r = new ReorderingDebugOutput(mDebugFormatter, debugOutput, 5 * mMaxGap)) {
        for (DiscordantReadSet drs : mWritten) {
          r.addRecord(drs);
        }
        for (DiscordantReadSet drs : mReadSets) {
          r.addRecord(drs);
        }
        // Also output a DRS corresponding to the current record only, useful for seeing how it overlaps the DRS it is merged into
        r.addRecord(new DiscordantReadSet(mTemplateName, mMaxGap, constraint));
      }
    }
    if (mParams.blockCompressed() && mParams.outputTabixIndex()) {
      indexTsv(dumpFile);
    }
  }

  static void processConstraint(final BreakpointConstraint constraint, final SortedSet<DiscordantReadSet> readSets, final String templateName, final int maxGap, SAMRecord record) {
    final Iterator<DiscordantReadSet> it = readSets.iterator();
    final List<DiscordantReadSet> overlap = new LinkedList<>();
    while (it.hasNext()) {
      final DiscordantReadSet drs = it.next();
      if (drs.belongs(constraint)) {
        it.remove();
        overlap.add(drs);
      }
    }
    final DiscordantReadSet newDrs;
    final int size = overlap.size();
    if (size == 0) {
      newDrs = new DiscordantReadSet(templateName, maxGap, constraint);
    } else if (size == 1) {
      newDrs = overlap.get(0);
      final boolean hasIntersection = newDrs.getIntersection() != null;
      newDrs.add(constraint);
      if (record != null && hasIntersection && newDrs.getIntersection() == null) {
        Diagnostic.developerLog("Intersection eliminated in " + newDrs + " with the addition of " + record.getSAMString());
      }
    } else {
      Diagnostic.developerLog("Merging multiple existing DRS: " + overlap);
      newDrs = new DiscordantReadSet(templateName, maxGap, constraint);
      for (final DiscordantReadSet over : overlap) {
        newDrs.addAll(over);
      }
    }
    if (record != null) {
      newDrs.add(record);
    }
    readSets.add(newDrs);

  }

  @Override
  protected void exec() throws IOException {
    try (OutputStream out = mParams.outStream(OUTPUT_FILENAME)) {
      mOut = out;
      try (OutputStream debugOutput = mParams.debugOutput() ? mParams.outStream(DEBUG_FILENAME) : null) {
        mDebugOutput = debugOutput;
        try {
          try (OutputStream bedOutput = mParams.bedOutput() ? mParams.outStream(BED_FILENAME) : null) {
            mBedOutput = bedOutput;
            super.exec();
            Diagnostic.userLog(mTotalDiscordantRecords + " discordant records.");
          }
        } finally {
          if (mDebugReordering != null) {
            mDebugReordering.close();
          }
        }
      }
    }
    if (mParams.blockCompressed() && mParams.outputTabixIndex()) {
      if (mDebugOutput != null) {
        indexTsv(mParams.outFile(DEBUG_FILENAME));
      }
      if (mBedOutput != null) {
        indexBed();
      }
      indexVcf();
    }
  }

  private void indexTsv(File file) throws IOException {
    try {
      new TabixIndexer(file).saveIndex(new DebugDiscordantOutputFormatter.DebugIndexerFactory());
    } catch (final UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + file + ": " + e.getMessage());
    }
  }

  private void indexBed() throws IOException {
    final Timer indexing = new Timer("BedIndex");
    indexing.start();
    final File file = mParams.outFile(BED_FILENAME);
    try {
      new TabixIndexer(file).saveBedIndex();
    } catch (final UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + file + ": " + e.getMessage());
    }
    indexing.stop();
    indexing.log();
  }

  private void indexVcf() throws IOException {
    final Timer indexing = new Timer("VcfIndex");
    indexing.start();
    final File file = mParams.outFile(OUTPUT_FILENAME);
    try {
      new TabixIndexer(file).saveVcfIndex();
    } catch (final UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + file + ": " + e.getMessage());
    }
    indexing.stop();
    indexing.log();
  }

  private class BreakpointPositionComparator implements Comparator<DiscordantReadSet> {
    @Override
    public int compare(DiscordantReadSet o1, DiscordantReadSet o2) {
      return new CompareHelper()
          .compare(readSetPosition(o1), readSetPosition(o2))
          .compare(o1.getUnion().getXLo(), o2.getUnion().getXLo())
          .compare(getConstraint(o1).getXHi(), getConstraint(o2).getXHi())
          .compare(getConstraint(o1).getYName(), getConstraint(o2).getYName())
          .compare(getConstraint(o1).getYLo(), getConstraint(o2).getYLo())
          .compare(getConstraint(o1).getYHi(), getConstraint(o2).getYHi())
          .result();
    }
  }

}
