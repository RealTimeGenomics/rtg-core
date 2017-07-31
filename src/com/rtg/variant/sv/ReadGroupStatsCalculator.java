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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.Environment;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Scans mapping results to calculate per-read group insert size
 * distribution and likelihood of a read being unmated/unmapped per
 * template position.
 *
 */
public final class ReadGroupStatsCalculator {

  static final String VERSION = "rgstats v1.3";

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
      mStats.computeIfAbsent(entry.getKey(), rgId -> new ReadGroupStats(rgId, entry.getValue().getTotalReference())).add(entry.getValue());
    }
  }


  /**
   * Set up the read group statistics map
   * @param header the SAM file header to create the map from
   */
  public void setupReadGroups(SAMFileHeader header) {
    final long refTotal = getRefTotal(header);
    for (final SAMReadGroupRecord rec : header.getReadGroups()) {
      mStats.computeIfAbsent(rec.getReadGroupId(), rgId -> new ReadGroupStats(rgId, refTotal));
    }
    if (mStats.size() == 0) {
      Diagnostic.warning("No read groups specified in SAM headers");
    }
  }

  private long getRefTotal(SAMFileHeader header) {
    long refTotal = 0;
    for (final SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
      refTotal += rec.getSequenceLength();
    }
    return refTotal;
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
  public void writeReadGroupStats(OutputStream out) throws IOException {
    long invalidCount = 0;
    long matedCount = 0;
    for (final ReadGroupStats stats : mStats.values()) {
      stats.calculate();
      Diagnostic.userLog(stats.toString());
      if (stats.isValid()) {
        break;
      } else {
        ++invalidCount;
        matedCount += stats.getMatedCount();
        Diagnostic.warning("Skipping read group with invalid values:" + StringUtils.LS + stats);
      }
    }
    if (invalidCount == mStats.size()) {
      if (matedCount == 0) {
        Diagnostic.warning("No properly mated alignments were available for read group statistics calculation");
      } else {
        Diagnostic.warning("No valid read group statistics could be calculated");
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

  /**
   * Add statistics from a single record.
   * @param record the SAM record
   */
  public void addRecord(SAMRecord record) {
    if (!record.getReadPairedFlag() || record.getReadUnmappedFlag() || record.getNotPrimaryAlignmentFlag() || record.isSecondaryOrSupplementary()) {
      return;
    }
    if (SamUtils.uniquelyMapped(record)) {
      final String rgId = ReadGroupUtils.getReadGroup(record);
      final ReadGroupStats stats = mStats.get(rgId);
      if (stats == null) {
        ++mNoReadGroupWarnings;
        if (mNoReadGroupWarnings <= 5) {
          Diagnostic.warning("Skipping record with no RG id");
        }
      } else {
        final Integer score = record.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES);
        if (score != null) {
          stats.addScore(score);
        }
        stats.addLength(record.getReadLength());
        if (record.getProperPairFlag()) {
          stats.addProper();
          final int tlen = Math.abs(record.getInferredInsertSize());
          final int gap = calculateGapSize(record);
          stats.addFragmentSize(tlen);
          stats.addGapSize(gap);
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
   * Add information from all records in the specified file
   * @param file the SAM file
   * @throws IOException if there was a problem reading the SAM data
   */
  public void addFile(File file) throws IOException {
    Diagnostic.userLog("Accumulating read group statistics for: " + file);
    mNoReadGroupWarnings = 0;
    try (RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(Collections.singletonList(file), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      setupReadGroups(it.header());
      while (it.hasNext()) {
        addRecord(it.next());
      }
      if (mNoReadGroupWarnings > 0) {
        Diagnostic.warning("Skipped " + mNoReadGroupWarnings + " records with no RG id in: " + file);
      }
    }
  }

  /**
   * Calculates gap size between SAM records
   * @param record SAM record
   * @return gap size
   */
  private static int calculateGapSize(SAMRecord record) {
    // XXX What does this do when the two reads overlap
    if (record.getMateAlignmentStart() > record.getAlignmentStart()) {
      return record.getMateAlignmentStart() - record.getAlignmentEnd();
    } else {
      return record.getAlignmentStart() - (record.getMateAlignmentStart() + record.getReadLength());
    }
  }

  void calculate(Collection<File> files, OutputStream out) throws IOException {
    for (final File f : files) {
      addFile(f);
    }
    writeReadGroupStats(out);
  }

}
