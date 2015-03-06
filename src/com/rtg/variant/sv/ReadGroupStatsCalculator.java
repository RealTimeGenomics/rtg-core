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
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SkipInvalidRecordsIterator;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * Scans mapping results to calculate per-read group insert size
 * distribution and likelihood of a read being unmated/unmapped per
 * template position.
 *
 */
public final class ReadGroupStatsCalculator {

  static final String VERSION = "rgstats v1.2";

  private static final String VERSION_HEADER = "#Version\t" + Environment.getVersion() + ", " + VERSION + StringUtils.LS;

  /**
   * Handles merging of multiple stats accumulation runs.
   */
  public static final class Merger {

    private final List<ReadGroupStatsCalculator> mCalculators = Collections.synchronizedList(new ArrayList<ReadGroupStatsCalculator>());
    /**
     * Creates an new ReadGroupStatsCalculator and adds it to the collective.
     * @return a ReadGroupStatsCalculator
     */
    public ReadGroupStatsCalculator createReadGroupStatsCalculator() {
      final ReadGroupStatsCalculator rgsc = new ReadGroupStatsCalculator();
      mCalculators.add(rgsc);
      return rgsc;
    }

    /**
     * Creates a ReadGroupStatsCalculator that has the combined statistics of all the other calculators.
     * @return a ReadGroupStatsCalculator
     */
    public ReadGroupStatsCalculator blend() {
      final ReadGroupStatsCalculator rgsc = new ReadGroupStatsCalculator();
      for (ReadGroupStatsCalculator irgsc : mCalculators) {
        rgsc.merge(irgsc);
      }
      rgsc.calculate();
      // TODO Clear mCalculators
      return rgsc;
    }

  }


  private int mNoReadGroupWarnings = 0;

  private final HashMap<String, ReadGroupStats> mStats;


  /**
   * Constructor
   */
  public ReadGroupStatsCalculator() {
    mStats = new HashMap<>();
  }


  /**
   * Merge contents of another read group statistics calculator into this one.
   * @param rgsc the read group statistics calculator to merge with this.
   */
  public void merge(ReadGroupStatsCalculator rgsc) {
    for (Entry<String, ReadGroupStats> entry : rgsc.mStats.entrySet()) {
      ReadGroupStats rgs = mStats.get(entry.getKey());
      if (rgs == null) {
        rgs = new ReadGroupStats(entry.getKey(), entry.getValue().getTotalReference());
        mStats.put(entry.getKey(), rgs);
      }
//      System.err.println("merge entry: " + entry.getKey() + " : " + entry.getValue().countsString());
//      System.err.println("merge: before: " + rgs.countsString());
      rgs.add(entry.getValue());
//      System.err.println("merge: after: " + rgs.countsString());
    }
  }


  /**
   * Set up the read group statistics map
   * @param header the SAM file header to create the map from
   */
  public void setupReadGroups(SAMFileHeader header) {
    long refTotal = 0;
    for (final SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
      refTotal += rec.getSequenceLength();
    }
    for (final SAMReadGroupRecord rec : header.getReadGroups()) {
      final String rgId = rec.getReadGroupId();
      final ReadGroupStats stats = mStats.get(rgId);
      if (stats == null) {
        mStats.put(rgId, new ReadGroupStats(rgId, refTotal));
      }
    }
    if (mStats.size() == 0) {
      throw new NoTalkbackSlimException("No read groups specified in SAM headers");
    }
  }

  /**
   * Returns the read group statistics for the given group
   * @param rgId the read group name
   * @return the stats for the desired group
   */
  public ReadGroupStats getStats(String rgId) {
    return mStats.get(rgId);
  }

  /**
   * Ensure all read group stats have been calculated
   */
  public void calculate() {
    for (final ReadGroupStats stats : mStats.values()) {
      stats.calculate();
      Diagnostic.userLog(stats.toString());
      if (stats.detectOverflow()) {
        throw new NoTalkbackSlimException("Overflow detected in read group statistics calculation for read group: " + stats.id());
      } else if (!stats.isValid()) {
        throw new NoTalkbackSlimException("Invalid read group statistics for read group: " + stats.id());
      }
    }
  }

  /**
   * Output calculated statistics
   * @param out output destination stream
   * @throws IOException when an IO exception occurs
   */
  public void dumpReadGroups(OutputStream out) throws IOException {
    long invalidCount = 0;
    long matedCount = 0;
    for (final ReadGroupStats stats : mStats.values()) {
      stats.calculate();
      if (stats.isValid()) {
        break;
      } else {
        invalidCount++;
        matedCount += stats.getMatedCount();
        Diagnostic.warning("Skipping read group with invalid values:" + StringUtils.LS + stats);
      }
    }
    if (invalidCount == mStats.size()) {
      if (matedCount == 0) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Read group statistics calculation requires properly mated data");
      } else {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "No valid read group statistics could be calculated");
      }
    }
    out.write(VERSION_HEADER.getBytes());
    out.write("#CL\t".getBytes());
    if (CommandLine.getCommandLine() != null) {
      out.write(CommandLine.getCommandLine().getBytes());
    } else {
      out.write("Internal".getBytes());
    }
    out.write((StringUtils.LS + ReadGroupStats.countsHeader() + StringUtils.LS).getBytes());
    for (final Map.Entry<String, ReadGroupStats> stats : mStats.entrySet()) {
      if (stats.getValue().isValid()) {
        out.write((stats.getValue().countsString() + StringUtils.LS).getBytes());
      }
    }
  }

  private void populateStats(File file, SamReader reader) throws IOException {
    try (RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(file.getPath(), reader)) {
      while (it.hasNext()) {
        final SAMRecord r = it.next();
        addRecord(r);
      }
      if (mNoReadGroupWarnings > 0) {
        Diagnostic.warning("Skipped " + mNoReadGroupWarnings + " records with no RG id");
      }
    }
  }

  /**
   * Add a single record.
   * @param record the SAM record
   */
  public void addRecord(SAMRecord record) {
    if (record.getReadUnmappedFlag() || !record.getReadPairedFlag()) {
      return;
    }
    Integer nh = SamUtils.getNHOrIH(record);
    if (nh == null) {
      //already should have warned in svprep
      nh = 1;
    }
    if (nh == 1) {
      final String rgId = ReadGroupUtils.getReadGroup(record);
      final ReadGroupStats stats = mStats.get(rgId);
      if (stats == null) {
        mNoReadGroupWarnings++;
        if (mNoReadGroupWarnings <= 5) {
          Diagnostic.warning("Skipping record with no RG id");
        }
      } else {
        stats.addScore(record.getIntegerAttribute("AS"));
        stats.addLength(record.getReadLength());
        if (record.getProperPairFlag()) {
          stats.addProper();
          stats.addFragmentSize(Math.abs(record.getInferredInsertSize()));
          stats.addGapSize(calculateGapSize(record));
        } else {
          if (record.getMateUnmappedFlag()) {
            stats.addUnmated();
          } else {
            stats.addDiscordant();
          }
        }
      }
    }
  }

  /**
   * Calculates gap size between SAM records
   * @param record SAM record
   * @return gap size
   */
  public static int calculateGapSize(SAMRecord record) {
    if (record.getMateAlignmentStart() > record.getAlignmentStart()) {
      return record.getMateAlignmentStart() - record.getAlignmentEnd();
    } else {
      return record.getAlignmentStart() - (record.getMateAlignmentStart() + record.getReadLength());
    }
  }

  /**
   * Performs the statistics calculation
   * @param files list of files to process
   * @param out output stream for standard output
   * @throws IOException when an IO error occurs
   */
  public void calculate(Collection<File> files, OutputStream out) throws IOException {
    for (final File f : files) {
      processFile(f);
    }
    dumpReadGroups(out);
  }

  private void processFile(File file) throws IOException {
    try (InputStream stream = FileUtils.createInputStream(file, false)) {
      try (SamReader reader = new SAMFileReader(stream)) {
        setupReadGroups(reader.getFileHeader());
        populateStats(file, reader);
      }
    } catch (final SAMException e) {
      if (e instanceof RuntimeIOException
          || e instanceof RuntimeEOFException) {
        throw new IOException(e.getMessage(), e);
      }
      throw new NoTalkbackSlimException(e, ErrorType.SAM_BAD_FORMAT_NO_FILE, e.getMessage());
    }
  }
}
