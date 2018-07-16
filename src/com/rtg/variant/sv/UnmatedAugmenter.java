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

package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.pairedend.InsertHelper;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceSequence;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.Pair;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;
import com.rtg.util.machine.PairOrientation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * Adds mating info to records of unique unmated hits with unique unmated hit for pair.
 * Ignores reads where one side does not have unique hit.
 */
public final class UnmatedAugmenter {

  /** Default filename for read group statistics */
  public static final String DEFAULT_RGSTATS_FILENAME = "rgstats.tsv";

  /**
   * Handles merging of multiple unmated augmenter runs.
   */
  public static final class Merger {

    private final List<UnmatedAugmenter> mAugmenters = Collections.synchronizedList(new ArrayList<UnmatedAugmenter>());
    private UnmatedAugmenter mBlended = null;

    /**
     * Creates an new UnmatedAugmenter and adds it to the collective.
     * @return a UnmatedAugmenter
     */
    public UnmatedAugmenter createUnmatedAugmenter() {
      final UnmatedAugmenter ua = new UnmatedAugmenter();
      mAugmenters.add(ua);
      return ua;
    }

    /**
     * Creates an UnmatedAugmenter that has the merged contents of the collective hash maps.
     * NOTE: This clears the collective as it merges them to allow garbage collection if no other references remain.
     * NOTE: Once this has been called, all future calls to it will return the same blended augmenter.
     * @return an UnmatedAugmenter
     */
    public synchronized UnmatedAugmenter blend() {
      if (mBlended == null) {
        mBlended = new UnmatedAugmenter();
        final Iterator<UnmatedAugmenter> iter = mAugmenters.iterator();
        while (iter.hasNext()) {
          final UnmatedAugmenter au = iter.next();
          mBlended.merge(au);
          au.mLeftSide.clear();
          au.mRightSide.clear();
          au.mRgMachineTypes.clear();
          iter.remove();
        }
      }
      return mBlended;
    }
  }

  // Contains minimal mate alignment information
  private static class CutRecord {
    final int mRefIndex;
    final int mAlignStart;
    final int mAlignEnd;
    final boolean mReverse;
    final int mAlignmentScore;

    CutRecord(SAMRecord r) {
      mRefIndex = r.getReferenceIndex();
      mAlignStart = r.getAlignmentStart() - 1;
      mAlignEnd = r.getAlignmentEnd();
      mReverse = r.getReadNegativeStrandFlag();
      mAlignmentScore = as(r);
    }

    private static int as(SAMRecord rec) {
      final Integer ii = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
      if (ii == null) {
        return Integer.MIN_VALUE;
      }
      return ii;
    }
  }

  private final HashMap<String, CutRecord> mLeftSide;
  private final HashMap<String, CutRecord> mRightSide;
  private final HashMap<String, MachineType> mRgMachineTypes;

  int mAugmentedUnmated = 0;
  int mAugmentedUnmapped = 0;

  /**
   * Constructor
   */
  public UnmatedAugmenter() {
    mLeftSide = new HashMap<>();
    mRightSide = new HashMap<>();
    mRgMachineTypes = new HashMap<>();
  }

  /**
   * Merge contents of another UnmatedAugmenter into this one.
   * @param other the other UnmatedAugmenter
   * @return a reference to this UnmatedAugmenter
   */
  private UnmatedAugmenter merge(UnmatedAugmenter other) {
    mLeftSide.putAll(other.mLeftSide);
    mRightSide.putAll(other.mRightSide);
    mRgMachineTypes.putAll(other.mRgMachineTypes);
    mAugmentedUnmapped += other.mAugmentedUnmapped;
    mAugmentedUnmated += other.mAugmentedUnmated;
    return this;
  }

  /**
   * Add record to internal hash maps used for augmenting if it is a unique unmated mapping without existing mate position information.
   * @param record SAM record to process.
   */
  public void addRecord(SAMRecord record) {
    if (isAugmentableUnmated(record)) {
      if (record.getFirstOfPairFlag()) {
        mLeftSide.put(record.getReadName(), new CutRecord(record));
      } else {
        mRightSide.put(record.getReadName(), new CutRecord(record));
      }
    }
  }

  private boolean isAugmentableUnmated(SAMRecord record) {
    if (!record.getReadPairedFlag() || record.getReadUnmappedFlag() || record.getProperPairFlag() || record.getMateReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
      return false;
    }
    final Integer nh = SamUtils.getNHOrIH(record);
    return nh != null && nh == 1;
  }

  /**
   * Adds machine type from a read group header line.
   * @param srgr read group record to get machine type from.
   */
  public void addMachineType(SAMReadGroupRecord srgr) {
    final MachineType mt = ReadGroupUtils.platformToMachineType(srgr, true);
    if (mt == null) {
      throw new NoTalkbackSlimException("Read group " + srgr.getId() + " does not contain a recognized platform");
    } else if (mt.orientation() == null) {
      throw new NoTalkbackSlimException("Platform " + srgr.getPlatform() + " is not supported (unknown expected read orientation)");
    }
    final String rgId = srgr.getReadGroupId();
    mRgMachineTypes.put(rgId, mt);
  }

  SAMFileWriter makeWriter(SAMFileHeader header, SamReader.Type type, OutputStream outStream, boolean presorted) {
    final SAMProgramRecord pg = new SAMProgramRecord(Constants.APPLICATION_NAME);
    pg.setProgramVersion(Environment.getVersion());
    if (CommandLine.getCommandLine() != null) {
      pg.setCommandLine(CommandLine.getCommandLine());
    } else {
      pg.setCommandLine("Internal");
    }
    SamUtils.addProgramRecord(header, pg);

    final SAMFileWriter writer;
    if (type == SamReader.Type.BAM_TYPE) {
      writer = new SAMFileWriterFactory().makeBAMWriter(header, presorted, outStream, true);
    } else {
      writer = new SAMFileWriterFactory().makeSAMWriter(header, presorted, outStream);
    }
    return writer;
  }

  /**
   * Augments a file containing only unmated alignment. Not reentrant.
   * @param file SAM/BAM file to augment
   * @param output SAM/BAM file to output.
   * @param calc statistics calculator to feed records through.
   * @throws IOException when an IO error occurs
   */
  void augmentMixed(File file, File output, ReadGroupStatsCalculator calc) throws IOException {
    // First pass, collect mated stats and augmentable unmated info
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(file), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      setupReadGroups(it.header(), calc);
      while (it.hasNext()) {
        final SAMRecord record = it.next();
        if (record.getProperPairFlag()) {
          calc.addRecord(record);
        }
        addRecord(record);
      }
    }

    calc.calculate(); // Ensure insert size stats have been computed

    // Second pass, augment unmated records (and collect their stats) and augment unmapped records (requires re-sorting)
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(file), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      try (SAMFileWriter writer = makeWriter(it.header(), SamUtils.getSamType(file), FileUtils.createOutputStream(output, FileUtils.isGzipFilename(file) || SamUtils.getSamType(file) == SamReader.Type.BAM_TYPE), false)) {
        while (it.hasNext()) {
          final SAMRecord record = it.next();
          if (!record.getProperPairFlag()) {
            if (record.getReadUnmappedFlag()) {
              updateUnmappedRecord(record, calc, null);
            } else {
              updateUnmatedRecord(record);
              calc.addRecord(record);
            }
          }
          writer.addAlignment(record);
        }
      }
    }
  }

  /**
   * Augments a file containing only unmated alignments. Not reentrant.
   * @param unmatedFile SAM/BAM file to augment
   * @param unmatedOutputFile SAM/BAM file to output.
   * @param calc statistics calculator to feed records through.
   * @throws IOException when an IO error occurs
   */
  public void augmentUnmated(File unmatedFile, File unmatedOutputFile, ReadGroupStatsCalculator calc) throws IOException {
    // First pass, collect augmentable unmated info
    try (RecordIterator<SAMRecord> it1 = new ThreadedMultifileIterator<>(Collections.singletonList(unmatedFile), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      setupReadGroups(it1.header(), calc);
      while (it1.hasNext()) {
        addRecord(it1.next());
      }
    }
    // Second pass, augment unmated records (and collect their stats)
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(unmatedFile), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      try (SAMFileWriter writer = makeWriter(it.header(), SamUtils.getSamType(unmatedFile), FileUtils.createOutputStream(unmatedOutputFile, FileUtils.isGzipFilename(unmatedFile) || SamUtils.getSamType(unmatedFile) == SamReader.Type.BAM_TYPE), true)) {
        while (it.hasNext()) {
          final SAMRecord r = it.next();
          updateUnmatedRecord(r);
          calc.addRecord(r);
          writer.addAlignment(r);
        }
      }
    }
  }

  private void setupReadGroups(SAMFileHeader header, ReadGroupStatsCalculator calc) {
    calc.setupReadGroups(header);
    for (final SAMReadGroupRecord srgr : header.getReadGroups()) {
      addMachineType(srgr);
    }
  }

  /**
   * Augments a file containing unmapped alignments. Not reentrant.
   * @param unmappedFile SAM/BAM file to augment
   * @param outputUnmappedFile SAM/BAM file to output.
   * @param calc statistics calculator to feed records through. <code>null</code> to ignore
   * @throws IOException when an IO error occurs
   */
  public void augmentUnmapped(File unmappedFile, File outputUnmappedFile, ReadGroupStatsCalculator calc) throws IOException {
    calc.calculate();
    // Augment the unmapped records (requires resorting)
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(unmappedFile), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      try (SAMFileWriter writer = makeWriter(it.header(), SamUtils.getSamType(unmappedFile), FileUtils.createOutputStream(outputUnmappedFile), false)) {
        while (it.hasNext()) {
          final SAMRecord r = it.next();
          if (r.getReadUnmappedFlag()) {
            updateUnmappedRecord(r, calc, null);
          }
          writer.addAlignment(r);
        }
      }
    }
  }

  /**
   * Update pair mapping information on an unmated SAM record based on the mapped location of it's mate..
   * @param record record to update.
   */
  public void updateUnmatedRecord(SAMRecord record) {
    if (isAugmentableUnmated(record)) {
      final CutRecord pair;
      if (record.getFirstOfPairFlag()) {
        pair = mRightSide.get(record.getReadName());
      } else {
        pair = mLeftSide.get(record.getReadName());
      }
      if (pair != null) { //NOTE: pair is only put in map if nh == 1
        mAugmentedUnmated++;
        record.setMateReferenceIndex(pair.mRefIndex);
        record.setMateAlignmentStart(pair.mAlignStart + 1);
        record.setMateNegativeStrandFlag(pair.mReverse);
        record.setMateUnmappedFlag(false);
        final int as = pair.mAlignmentScore;
        if (as >= 0) {
          record.setAttribute(SamUtils.ATTRIBUTE_MATE_ALIGNMENT_SCORE, as);
        }
        record.setAttribute(SamUtils.ATTRIBUTE_MATE_END, pair.mAlignEnd);
        if (pair.mRefIndex == record.getReferenceIndex()) {
          final int tlen = InsertHelper.tlen(record.getFirstOfPairFlag(), record.getAlignmentStart(), record.getAlignmentEnd() - record.getAlignmentStart() + 1, pair.mAlignStart + 1, pair.mAlignEnd - pair.mAlignStart);
          record.setInferredInsertSize(tlen);
        } else {
          record.setInferredInsertSize(0);
        }
      }
    }
  }
  /**
   * Updates an unmapped record with pair mapping information.
   * @param record SAM record to update
   * @param calc a read group statistics calculator upon which calc.calculate() has been called.
   * @param referenceGenome reference genome for reference, or null if unknown or PAR region support is not required
   */
  public void updateUnmappedRecord(SAMRecord record, ReadGroupStatsCalculator calc, ReferenceGenome referenceGenome) {
    assert record.getReadUnmappedFlag();
    final CutRecord mate;
    if (record.getFirstOfPairFlag()) {
      mate = mRightSide.get(record.getReadName());
    } else {
      mate = mLeftSide.get(record.getReadName());
    }
    if (mate != null) { //other side was mapped NOTE: mate only in map if nh==1
      mAugmentedUnmapped++;
      record.setMateReferenceIndex(mate.mRefIndex);
      record.setMateAlignmentStart(mate.mAlignStart + 1);
      record.setMateNegativeStrandFlag(mate.mReverse);
      record.setMateUnmappedFlag(false);
      record.setAttribute(SamUtils.ATTRIBUTE_MATE_END, mate.mAlignEnd);

      //make up some start positions etc for this read.
      final String rg = ReadGroupUtils.getReadGroup(record);
      final MachineType mt = mRgMachineTypes.get(rg);

      if (mt != null) {
        final MachineOrientation machineOrientation = mt.orientation();
        final PairOrientation mateOrientation;
        if (mate.mReverse) {
          mateOrientation = record.getFirstOfPairFlag() ? PairOrientation.R2 : PairOrientation.R1;
        } else {
          mateOrientation = record.getFirstOfPairFlag() ? PairOrientation.F2 : PairOrientation.F1;
        }
        final PairOrientation readOrientation = machineOrientation.getMateOrientation(mateOrientation);
        if (readOrientation != null) {
          record.setReferenceIndex(mate.mRefIndex);
          if (PairOrientation.F1 == readOrientation || PairOrientation.F2 == readOrientation) {
            record.setReadNegativeStrandFlag(false);
          } else {
            record.setReadNegativeStrandFlag(true);
            //need to RC the read and rev the qualities.

            final byte[] b = new byte[record.getReadBases().length];
            DnaUtils.reverseComplement(record.getReadBases(), b, b.length);
            record.setReadBases(b);

            final byte[] q = new byte[record.getBaseQualities().length];
            for (int i = 0; i < q.length; ++i) {
              q[q.length - 1 - i] = record.getBaseQualities()[i];
            }
            record.setBaseQualities(q);
          }

          final ReadGroupStats rgstats = calc.getStats(rg);
          final int thisFragmentLength;
          if (rgstats != null) {
            thisFragmentLength = (int) rgstats.fragmentMean();
          } else if (record.getReadGroup().getPredictedMedianInsertSize() != null) {
            thisFragmentLength = record.getReadGroup().getPredictedMedianInsertSize();
          } else {
            throw new NoTalkbackSlimException("Could not determine mean fragment size for placing unmapped records.");
          }

          setPlacement(record, referenceGenome, machineOrientation, mateOrientation, mate.mAlignStart, mate.mAlignEnd, thisFragmentLength);
        }
      }
    }
  }

  static void setPlacement(SAMRecord record, ReferenceGenome referenceGenome, MachineOrientation machineOrientation, PairOrientation mateOrientation, int mateStart, int mateEnd, int fragmentLength) {
    int alignmentStart;
    if (machineOrientation.isMateUpstream(mateOrientation)) {
      alignmentStart = Math.max(0, mateStart + fragmentLength - record.getReadLength());
    } else {
      alignmentStart = Math.max(0, mateEnd - fragmentLength);
    }
    final int refLength = record.getHeader().getSequence(record.getReferenceIndex()).getSequenceLength();
    if (referenceGenome != null) {
      // If our projected position is in the duplicated PAR region, port to the canonical location
      final ReferenceSequence rs = referenceGenome.sequence(record.getReferenceName());
      if (rs.hasDuplicates()) {
        for (Pair<RegionRestriction, RegionRestriction> dup : rs.duplicates()) {
          if (dup.getB().contains(record.getReferenceName(), Math.min(refLength, alignmentStart))) {
            alignmentStart = dup.getA().getStart() + (alignmentStart - dup.getB().getStart());
            record.setReferenceName(dup.getA().getSequenceName());
          }
        }
      }
    }
    alignmentStart = Math.min(refLength, alignmentStart);
    record.setAlignmentStart(alignmentStart + 1);
  }

}
