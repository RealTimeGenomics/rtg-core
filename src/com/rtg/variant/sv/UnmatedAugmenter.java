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
import com.rtg.sam.SamUtils;
import com.rtg.sam.SkipInvalidRecordsIterator;
import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineType;
import com.rtg.util.machine.PairOrientation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

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


  private final HashMap<String, CutRecord> mLeftSide;
  private final HashMap<String, CutRecord> mRightSide;
  private final HashMap<String, MachineType> mRgMachineTypes;

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
    return this;
  }

  private void populateMaps(File file) throws IOException {
    try (SAMFileReader reader = new SAMFileReader(FileUtils.createFileInputStream(file, false))) {
      try (RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(file.getPath(), reader)) {
        while (it.hasNext()) {
          processRecord(it.next());
        }
      }
    }
  }

  /**
   * Process the given SAM record, will be added to internal hash maps if it is a unique unmated mapping.
   * @param record SAM record to process.
   */
  public void processRecord(SAMRecord record) {
    final Integer nh = SamUtils.getNHOrIH(record);
    if (record.getReadPairedFlag() && !record.getProperPairFlag() && nh != null && nh == 1) {
      if (record.getFirstOfPairFlag()) {
        mLeftSide.put(record.getReadName(), new CutRecord(record));
      } else {
        mRightSide.put(record.getReadName(), new CutRecord(record));
      }
    }
  }

  private static int as(SAMRecord rec) {
    final Integer ii = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
    if (ii == null) {
      return Integer.MIN_VALUE;
    }
    return ii;
  }

  private static class CutRecord {
    final int mRefIndex;
    final int mAlignStart;
    final int mAlignEnd;
    final boolean mReverse;
    final int mAlignmentScore;

    CutRecord(SAMRecord r) {
      mRefIndex = r.getReferenceIndex();
      mAlignStart = r.getAlignmentStart();
      mAlignEnd = r.getAlignmentEnd();
      mReverse = r.getReadNegativeStrandFlag();
      mAlignmentScore = as(r);
    }
  }

  private SAMFileHeader constructHeader(SAMFileReader reader, boolean processMachineTypes) {
    final SAMFileHeader header = reader.getFileHeader();
    final SAMProgramRecord pg = new SAMProgramRecord(Constants.APPLICATION_NAME);
    pg.setProgramVersion(Environment.getVersion());
    if (CommandLine.getCommandLine() != null) {
      pg.setCommandLine(CommandLine.getCommandLine());
    } else {
      pg.setCommandLine("Internal");
    }
    SamUtils.addProgramRecord(header, pg);

    if (processMachineTypes) {
      for (final SAMReadGroupRecord srgr : header.getReadGroups()) {
        addMachineType(srgr);
      }
    }
    return header;
  }

  /**
   * Adds machine type from a read group header line.
   * @param srgr read group record to get machine type from.
   */
  public void addMachineType(SAMReadGroupRecord srgr) {
    final MachineType mt = ReadGroupUtils.platformToMachineType(srgr, true);
    if (mt == null) {
      throw new NoTalkbackSlimException("Read group with platform specified required");
    } else if (mt.orientation() == null) {
      throw new NoTalkbackSlimException("Unsupported platform: " + srgr.getPlatform());
    }
    final String rgId = srgr.getReadGroupId();
    mRgMachineTypes.put(rgId, mt);
  }

  private SAMFileWriter makeWriter(SAMFileReader reader, OutputStream outStream, boolean processMachineTypes, boolean presorted) {
    final SAMFileHeader header = constructHeader(reader, processMachineTypes);
    final SAMFileWriter writer;
    if (reader.isBinary()) {
      writer = new SAMFileWriterFactory().makeBAMWriter(header, presorted, outStream, true);
    } else {
      writer = new SAMFileWriterFactory().makeSAMWriter(header, presorted, outStream);
    }
    return writer;
  }

  /**
   * Performs the augmentation. Not reentrant.
   * @param unmappedFile SAM/BAM file to augment
   * @param outputUnmappedFile SAM/BAM file to output.
   * @param zipInput whether input is SAM and <code>gzipped</code>
   * @param calc statistics calculator to feed records through. <code>null</code> to ignore
   * @throws IOException when an IO error occurs
   */
  public void augmentUnmapped(File unmappedFile, File outputUnmappedFile, boolean zipInput, ReadGroupStatsCalculator calc) throws IOException {
    calc.calculate();
    try (SAMFileReader reader = new SAMFileReader(FileUtils.createFileInputStream(unmappedFile, false))) {
      try (OutputStream outStream = FileUtils.createOutputStream(outputUnmappedFile, zipInput)) {
        try (SAMFileWriter writer = makeWriter(reader, outStream, true, false)) {
          final RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(unmappedFile.getPath(), reader, true);
          while (it.hasNext()) {
            final SAMRecord r = it.next();

            updateUnmappedRecord(r, calc, null);

            writer.addAlignment(r);
          }
        }
      }
    }
  }

  /**
   * Updates an unmapped record with pair mapping information.
   * Assumes calc.calculate() has been called prior to this message.
   * @param record SAM record to update
   * @param calc a read group statistics calculator
   * @param referenceGenome reference genome for reference, or null if unknown
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
      record.setMateReferenceIndex(mate.mRefIndex);
      record.setMateAlignmentStart(mate.mAlignStart);
      record.setMateNegativeStrandFlag(mate.mReverse);
      record.setMateUnmappedFlag(false);
      record.setAttribute(SamUtils.ATTRIBUTE_MATE_END, mate.mAlignEnd);

      //make up some start positions etc for this read.
      final String rg = ReadGroupUtils.getReadGroup(record);
      final MachineType mt = mRgMachineTypes.get(rg);

      if (mt != null) {
        final PairOrientation mateOrientation;
        if (mate.mReverse) {
          mateOrientation = record.getFirstOfPairFlag() ? PairOrientation.R2 : PairOrientation.R1;
        } else {
          mateOrientation = record.getFirstOfPairFlag() ? PairOrientation.F2 : PairOrientation.F1;
        }
        final PairOrientation po = mt.orientation().getMateOrientation(mateOrientation);
        if (po != null) {
          record.setReferenceIndex(mate.mRefIndex);
          if (PairOrientation.F1.equals(po) || PairOrientation.F2.equals(po)) {
            record.setReadNegativeStrandFlag(false);
          } else {
            record.setReadNegativeStrandFlag(true);
            //need to RC the read and rev the qualities.

            final byte[] b = new byte[record.getReadBases().length];
            DnaUtils.reverseComplement(record.getReadBases(), b, b.length);
            record.setReadBases(b);

            final byte[] q = new byte[record.getBaseQualities().length];
            for (int i = 0; i < q.length; i++) {
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
          int alignmentStart;
          if (mt.orientation().isMateUpstream(mateOrientation)) {
            alignmentStart = Math.max(1, mate.mAlignStart + thisFragmentLength - record.getReadLength());
          } else {
            alignmentStart = Math.max(1, mate.mAlignEnd - 1 - thisFragmentLength);
          }
          final int refLength = record.getHeader().getSequence(record.getReferenceIndex()).getSequenceLength();
          if (referenceGenome != null) {
            //get ReferenceSequence and check if we are in the duplicated par region
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
          record.setAlignmentStart(Math.min(refLength, alignmentStart));
        }
      }
    }
  }

  /**
   * Updates read group statistics with information from a mated SAM file.
   * @param matedFile mated SAM or BAM file
   * @param calc ReadGroupStatsCalculator object to update
   * @throws IOException in an error occurs while reading file
   */
  public static void populateReadGroupStats(File matedFile, ReadGroupStatsCalculator calc) throws IOException {
    try (SAMFileReader reader = new SAMFileReader(FileUtils.createFileInputStream(matedFile, false))) {
      calc.setupReadGroups(reader.getFileHeader());

      //assume any warnings from input are reported the first time we processed so set Skipping iterator to silent mode.
      final RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(matedFile.getPath(), reader, true);
      while (it.hasNext()) {
        final SAMRecord r = it.next();
        calc.addRecord(r);
      }
    }
  }

  /**
   * Performs the augmentation. Not reentrant.
   * @param unmatedFile SAM/BAM file to augment
   * @param unmatedOutputFile SAM/BAM file to output.
   * @param zipInput whether input is SAM and <code>gzipped</code>
   * @param calc statistics calculator to feed records through. <code>null</code> to ignore
   * @throws IOException when an IO error occurs
   */
  public void augmentUnmated(File unmatedFile, File unmatedOutputFile, boolean zipInput, ReadGroupStatsCalculator calc) throws IOException {
    populateMaps(unmatedFile);
    try (SAMFileReader reader = new SAMFileReader(FileUtils.createFileInputStream(unmatedFile, false))) {
      if (calc != null) {
        calc.setupReadGroups(reader.getFileHeader());
      }

      try (SAMFileWriter writer = makeWriter(reader, FileUtils.createOutputStream(unmatedOutputFile, zipInput || reader.isBinary()), false, true)) {

        //assume any warnings from input are reported the first time we processed so set Skipping iterator to silent mode.
        final RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(unmatedFile.getPath(), reader, true);
        while (it.hasNext()) {
          final SAMRecord r = it.next();

          updateUnmatedRecord(r);

          if (calc != null) {
            calc.addRecord(r);
          }
          writer.addAlignment(r);
        }
      }
    }
  }

  /**
   * Update pair mapping information on an unmated SAM record.
   * @param record record to update.
   */
  public void updateUnmatedRecord(SAMRecord record) {
    if (!record.getReadUnmappedFlag()) {
      final Integer nh = SamUtils.getNHOrIH(record);
      if (nh == null) {
        throw new NoTalkbackSlimException("No NH or IH value in mapped record: " + record);
      }
      if (nh == 1 && record.getReadPairedFlag() && !record.getProperPairFlag()) {
        final CutRecord pair;
        if (record.getFirstOfPairFlag()) {
          pair = mRightSide.get(record.getReadName());
        } else {
          pair = mLeftSide.get(record.getReadName());
        }
        if (pair != null) { //NOTE: pair is only put in map if nh == 1
          record.setMateReferenceIndex(pair.mRefIndex);
          record.setMateAlignmentStart(pair.mAlignStart);
          record.setMateNegativeStrandFlag(pair.mReverse);
          record.setMateUnmappedFlag(false);
          final int as = pair.mAlignmentScore;
          if (as >= 0) {
            record.setAttribute(SamUtils.ATTRIBUTE_MATE_ALIGNMENT_SCORE, as);
          }
          record.setAttribute(SamUtils.ATTRIBUTE_MATE_END, pair.mAlignEnd);
          if (pair.mRefIndex == record.getReferenceIndex()) {
            final int tlen = InsertHelper.tlen(record.getFirstOfPairFlag(), record.getAlignmentStart(), record.getAlignmentEnd() - record.getAlignmentStart() + 1, pair.mAlignStart, pair.mAlignEnd - pair.mAlignStart + 1);
            record.setInferredInsertSize(tlen);
          } else {
            record.setInferredInsertSize(0);
          }
        }
      }
    }
  }
}
