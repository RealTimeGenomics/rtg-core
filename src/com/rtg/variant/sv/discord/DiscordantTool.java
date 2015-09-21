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
package com.rtg.variant.sv.discord;


import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.sam.BamIndexer;
import com.rtg.sam.SamIteratorTask;
import com.rtg.sam.SamUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.CompareHelper;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.sv.ReadGroupStats;
import com.rtg.vcf.VcfRecord;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 */
public class DiscordantTool extends SamIteratorTask<DiscordantToolParams, DiscordantToolStatistics> {
  private static final String DISCORD_RECORDS_BAM = "discord.records.bam";
  static final String INDIVIDUAL = "individual";
  static final String NONE = "none";
  static final String DUMP_RECORDS = NONE; //System.getProperty("discord.dump.records", NONE);

  static final String OUTPUT_FILENAME = "discordant_pairs.vcf";
  static final String DEBUG_FILENAME = "debug_discordant_pairs.txt";
  static final String BED_FILENAME = "discordant_pairs.bed";
  static final String DT_OUTPUT_VERSION = "1";


  private int mMaxGap;
  protected final SortedSet<DiscordantReadSet> mReadSets = new TreeSet<>(new DiscordantReadSet.FlushPositionComparator());

  private String mSampleName;
  final Map<String, MachineOrientation> mMachineOrientations = new HashMap<>();

  private final VcfDiscordantOutputFormatter mFormatter;
  private final DebugDiscordantOutputFormatter mDebugFormatter;
  private final BedDiscordantOutputFormatter mBedFormatter;
  private ReorderingDebugOutput mDebugReordering;
  private final TreeSet<DiscordantReadSet> mOutputBuffer = new TreeSet<>(new BreakpointPositionComparator());

  private OutputStream mOut = null;
  private OutputStream mDebugOutput = null;
  private OutputStream mBedOutput = null;
  private SmartBedWriter mBedWriter = null;
  private SmartVcfWriter mVcfWriter = null;

  private SAMFileHeader mSamHeader = null;
  private long mTotalDiscordantRecords = 0;
  private SAMFileWriter mBamDumpWriter = null; //used for system property discord.dump.records
  private File mBamDump = null;

  protected DiscordantTool(DiscordantToolParams params, OutputStream defaultOutput) throws IOException {
    super(params, defaultOutput, new DiscordantToolStatistics(params.directory()), params.filterParams());
    mFormatter = new VcfDiscordantOutputFormatter(mGenomeSequences);
    mBedFormatter = params.bedOutput() ? new BedDiscordantOutputFormatter() : null;
    mDebugFormatter =  params.debugOutput() ? new DebugDiscordantOutputFormatter() : null;
  }

  BreakpointConstraint getConstraint(DiscordantReadSet rs) {
    return rs.getIntersection() == null ? rs.getUnion() : rs.getIntersection();
  }
  int readSetPosition(DiscordantReadSet rs) {
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
    Diagnostic.developerLog(DEBUG_FILENAME);
    if (mBamDumpWriter != null) {
      mBamDumpWriter.close();
      try {
        BamIndexer.saveBamIndex(mBamDump, new File(mBamDump.getParent(), mBamDump.getName() + BamIndexer.BAM_INDEX_EXTENSION));
      } catch (final UnindexableDataException e) {
        Diagnostic.warning("Cannot produce index for: " + mBamDump + ": " + e.getMessage());
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

  void flush(DiscordantReadSet drs) {
    if (drs.getCounts() >= mParams.minBreakpointDepth() && (!mParams.intersectionOnly() || drs.getIntersection() != null)) {
      mOutputBuffer.add(drs);
    }
  }
  void output() throws IOException {
    while (!mOutputBuffer.isEmpty()) {
      final DiscordantReadSet first = mOutputBuffer.pollFirst();
      writeReadSet(first, -1, -1);
    }
  }

  void writeReadSet(DiscordantReadSet drs, int coverage, double ambiguous) throws IOException {
    if (mParams.debugOutput() || !DUMP_RECORDS.equals(NONE)) {
      final String debugRep = mDebugFormatter.format(drs);
      if (mParams.debugOutput()) {
        mDebugReordering.addRecord(drs);
      }
      //mDebugOutput.write(StringUtils.LS.getBytes());
      if (!DUMP_RECORDS.equals(NONE)) {
        final SAMFileHeader header = mSamHeader.clone();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addComment(debugRep.trim());
        final SAMFileWriter sw;
        final File bamoutput;
        if (DUMP_RECORDS.equals(INDIVIDUAL)) {
          bamoutput = mParams.file(drs.getSequenceName() + ":" + drs.unionPosition() + ".bam");
          sw = new SAMFileWriterFactory().makeBAMWriter(header, false, FileUtils.createOutputStream(bamoutput, false, false), false);
        } else {
          bamoutput = null;
          if (mBamDumpWriter == null) {
            mBamDump = mParams.file(DISCORD_RECORDS_BAM);
            mBamDumpWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, FileUtils.createOutputStream(mBamDump, false, false), false);
          }
          sw = mBamDumpWriter;
        }

        // Output to uncompressed SAM so we can load into IGV without indexing?
        for (final SAMRecord rec : drs.getRecords()) {
          sw.addAlignment(rec);
        }
        if (DUMP_RECORDS.equals(INDIVIDUAL)) {
          sw.close();
          try {
            BamIndexer.saveBamIndex(bamoutput, new File(bamoutput.getParent(), bamoutput.getName() + BamIndexer.BAM_INDEX_EXTENSION));
          } catch (final UnindexableDataException e) {
            Diagnostic.warning("Cannot produce index for: " + bamoutput + ": " + e.getMessage());
          }
        }
      }
    }
    final VcfRecord rec = mFormatter.vcfRecord(drs, coverage, ambiguous);
    mVcfWriter.addRecord(rec);
    mStatistics.tallyDiscordantReadSet(drs);
    if (mBedOutput != null) {
      mBedWriter.addRecord(mBedFormatter.format(drs));
    }
  }

  @Override
  public void init(SAMFileHeader header) throws IOException {
    if (header.getReadGroups() == null || header.getReadGroups().size() == 0) {
      throw new NoTalkbackSlimException("Read group header is not present in the SAM file");
    }
    mSamHeader = header;
    int maxGap = 0;
    boolean first = true;
    for (final SAMReadGroupRecord srgr : header.getReadGroups()) {
      if (srgr.getPlatform() == null) {
        throw new NoTalkbackSlimException("Platform in Read Group header " + srgr.getReadGroupId() + " was not present");
      }
      if (first) {
        first = false;
        mSampleName = srgr.getSample();
        assert mSampleName != null;
      }
      if (mSampleName != null && !mSampleName.equals(srgr.getSample())) {
        throw new NoTalkbackSlimException("Sample name does not match " + mSampleName + ", " + srgr.getSample());
      }
      final ReadGroupStats rgs = mParams.readGroupStatistics().get(srgr.getId());
      if (rgs == null) {
        throw new NoTalkbackSlimException("Read Group " + srgr.getId() + " did not have corresponding rgstats data supplied");
      }
      final int gap = BreakpointConstraint.gapMax(rgs);
      maxGap = Math.max(maxGap, gap);
      final String platform = srgr.getPlatform().toUpperCase(Locale.ROOT);
      if (MachineType.ILLUMINA_PE.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.ILLUMINA_PE.orientation());
      } else if (MachineType.COMPLETE_GENOMICS.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.COMPLETE_GENOMICS.orientation());
      } else if (MachineType.COMPLETE_GENOMICS_2.compatiblePlatform(platform)) {
        mMachineOrientations.put(srgr.getId(), MachineType.COMPLETE_GENOMICS_2.orientation());
      } else {
        throw new NoTalkbackSlimException("Unsupported Platform: " + srgr.getPlatform());
      }
    }
    mMaxGap = maxGap;
    if (mDebugOutput != null) {
      mDebugOutput.write(mDebugFormatter.header().getBytes());
      mDebugReordering = new ReorderingDebugOutput(mDebugFormatter, mDebugOutput, 5 * mMaxGap);
    }
    if (mBedOutput != null) {
      mBedOutput.write(mBedFormatter.header().getBytes()); // Maybe move this into the bedwriter
      mBedWriter = new SmartBedWriter(mBedOutput);
    }
    mVcfWriter = new SmartVcfWriter(mOut, mFormatter.header(header, mSampleName, false, false));
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
    final Integer nh = SamUtils.getNHOrIH(rec);
    if (nh == null || nh < 2)  {
      final String id = rec.getReadGroup().getId();
      final MachineOrientation mo = mMachineOrientations.get(id);
      final ReadGroupStats rgs = mParams.readGroupStatistics().get(id);
      final BreakpointConstraint constraint = new BreakpointConstraint(rec, mo, rgs);
      if (constraint.isConcordant()) {
        return true;
      }
      processConstraint(constraint, mReadSets, mTemplateName, mMaxGap, rec);
      mTotalDiscordantRecords++;
    }
    flush(0, rec.getAlignmentStart());
    return true;
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
      newDrs.add(constraint);
    } else {
      newDrs = new DiscordantReadSet(templateName, maxGap, constraint);
      for (final DiscordantReadSet over : overlap) {
        newDrs.addAll(over);
      }
    }
    if (!DUMP_RECORDS.equals(NONE)) {
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
          if (mDebugOutput != null) {
            mDebugReordering.close();
          }
        }
      }
    }
    if (mParams.blockCompressed() && mParams.outputTabixIndex()) {
      if (mBedOutput != null) {
        final Timer indexing = new Timer("BedIndex");
        indexing.start();
        final File file = mParams.outFile(BED_FILENAME);
        try {
          new TabixIndexer(file, new File(file.getParent(), file.getName() + TabixIndexer.TABIX_EXTENSION)).saveBedIndex();
        } catch (final UnindexableDataException e) {
          Diagnostic.warning("Cannot produce TABIX index for: " + file + ": " + e.getMessage());
        }
        indexing.stop();
        indexing.log();
      }
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
  }

  private class BreakpointPositionComparator implements Comparator<DiscordantReadSet> {
    @Override
    public int compare(DiscordantReadSet o1, DiscordantReadSet o2) {
      return new CompareHelper()
          .compare(readSetPosition(o1), readSetPosition(o2))
          .compare(o1.getUnion().getX(), o2.getUnion().getX())
          .compare(getConstraint(o1).getZ(), getConstraint(o2).getZ())
          .compare(getConstraint(o1).getYName(), getConstraint(o2).getYName())
          .compare(getConstraint(o1).getY(), getConstraint(o2).getY())
          .compare(getConstraint(o1).getW(), getConstraint(o2).getW())
          .result();
    }
  }

}
