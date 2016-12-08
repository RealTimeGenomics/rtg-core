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


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.NormalDistribution;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Encapsulates priors and distributions for alignments within a read group.
 *
 */
public class ReadGroupStats {

  private static final double TINY = 0.0000001;

  /**
   * Reads a set of read group name remappings from source name to
   * destination name. This is used to merge multiple read groups into
   * one. Format is tab-separated lines consisting of input read group
   * name, output read group name.
   *
   * @param relabelFile the file containing the group name mappings.
   * @return a map from input read group name to output read group name.
   * @exception IOException if an error occurs.
   */
  static Map<String, String> loadRelabelFile(File relabelFile) throws IOException {
    final HashMap<String, String> remap = new HashMap<>();
    Diagnostic.userLog("Loading read group relabelling file");
    try (BufferedReader br = new BufferedReader(new FileReader(relabelFile))) {
      String line;
      while ((line = br.readLine()) != null) {
        if (!line.startsWith("#") && line.length() > 0) {
          final String[] parts = line.split("\t");
          if (parts.length != 2) {
            throw new IOException("Expected input_name<tab>output_name on line: " + line);
          }
          remap.put(parts[0], parts[1]);
        }
      }
      if (remap.size() == 0) {
        throw new IOException("No read group relabellings contained in file " + relabelFile);
      }
    }
    return remap;
  }

  /**
   * Reads a set of read group priors.
   *
   * @param remap optional read group id remapping (used to merge read group ids)
   * @param files the files containing the priors.
   * @return a map from read group name to a ReadGroupStats object.
   * @exception IOException if an error occurs.
   */
  static Map<String, ReadGroupStats> loadReadGroupStats(Map<String, String> remap, File... files) throws IOException {
    final HashMap<String, ReadGroupStats> rgStats = new HashMap<>();
    Diagnostic.userLog("Loading read group statistics");
    for (final File statsFile : files) {
      try (BufferedReader br = new BufferedReader(new FileReader(statsFile))) {
        String line;
        while ((line = br.readLine()) != null) {
          if (line.length() > 0) {
            if (line.startsWith("#")) {
              if (line.startsWith("#Version") && !line.contains(ReadGroupStatsCalculator.VERSION)) {
                throw new NoTalkbackSlimException("Unsupported rgstats version: " + line + " - rerun svprep");
              }
            } else {
              try {
                final ReadGroupStats stats = new ReadGroupStats(remap, line);
                if (rgStats.containsKey(stats.id())) {
                  rgStats.get(stats.id()).add(stats);
                } else {
                  rgStats.put(stats.id(), stats);
                }
              } catch (final IllegalArgumentException e) {
                throw new IOException(e.getMessage());
              }
            }
          }
        }
        if (rgStats.size() == 0) {
          throw new IOException("No read group statistics contained in file " + statsFile);
        }
      }
    }
    for (final ReadGroupStats stats : rgStats.values()) {
      stats.calculate();
      Diagnostic.userLog(stats.toString());
      if (stats.detectOverflow()) {
        throw new NoTalkbackSlimException("Overflow detected in read group statistics calculation for read group: " + stats.id());
      } else if (!stats.isValid()) {
        throw new NoTalkbackSlimException("Invalid read group statistics for read group: " + stats.id());
      }
    }
    return rgStats;
  }


  private String mReadGroup;
  private int mMaxAlignment;
  private long mTotalReference;
  private NormalDistribution mLengthDist = new NormalDistribution();
  private NormalDistribution mFragmentSizeDist = new NormalDistribution();
  private NormalDistribution mGapSizeDist = new NormalDistribution();
  private long mProper = 0;
  private long mDiscordant = 0;
  private long mUnmated = 0;

  // XXX Assume 8 percent of the read length may have been erroneously aligned across a breakpoint due to alignment penalties.
  // This should really be passed in somewhere. We used to assume that mAlignmentStartIgnored == mMaxAlignment, but since our
  // alignment penalties have changed, this is no longer true.
  private static final int ALIGNMENT_START_IGNORED_FRACTION = 8;

  // Derived from counts stored above
  private double mMeanLength;
  private double mFragmentMean;
  private double mFragmentStdDev;
  private double mGapMean;
  private double mGapStdDev;
  private double mProperRate;
  private double mProperRandomRate;
  private double mDiscordantRate;
  private double mUnmatedRate;
  private int mAlignmentStartIgnored;

  protected ReadGroupStats(String id, long totalReference) {
    mReadGroup = id;
    mTotalReference = totalReference;
  }

  /**
   * Constructs a ReadGroupStats by parsing the fields out of a
   * tab-separated line
   * @param remap optional read group id remapping (used to merge read group ids)
   * @param statsLine the line to parse
   * @throws IllegalArgumentException if the line is malformed
   */
  ReadGroupStats(Map<String, String> remap, String statsLine) {
    final String[] parts = statsLine.split("\t");
    if (parts.length != 15) {
      throw new IllegalArgumentException("Incorrect number of fields on line: " + statsLine);
    }
    int fnum = 0;
    mReadGroup = parts[fnum++];
    if ((remap != null) && remap.containsKey(mReadGroup)) { // If no mapping provided, keep original id
      mReadGroup = remap.get(mReadGroup);
    }
    try {
      mTotalReference = Long.parseLong(parts[fnum++]);
      mLengthDist = new NormalDistribution(Integer.parseInt(parts[fnum++]), Double.parseDouble(parts[fnum++]), Double.parseDouble(parts[fnum++]));
      mFragmentSizeDist = new NormalDistribution(Integer.parseInt(parts[fnum++]), Double.parseDouble(parts[fnum++]), Double.parseDouble(parts[fnum++]));
      mGapSizeDist = new NormalDistribution(Integer.parseInt(parts[fnum++]), Double.parseDouble(parts[fnum++]), Double.parseDouble(parts[fnum++]));
      mMaxAlignment = Integer.parseInt(parts[fnum++]);
      mProper = Long.parseLong(parts[fnum++]);
      mDiscordant = Long.parseLong(parts[fnum++]);
      mUnmated = Long.parseLong(parts[fnum++]);
      calculate();
    } catch (final NumberFormatException e) {
      throw new IllegalArgumentException("Badly formatted field " + fnum + " on line: " + statsLine);
    }
  }

  /**
   * Only to be used by tests.
   *
   * @param id identifier.
   * @param meanLength mean length.
   * @param fragmentMean mean of fragment length
   * @param fragmentStdDev standard deviation of fragment length
   * @param gapMean mean of gap length.
   * @param gapStdDev standard deviation of gap length.
   * @param maxAlignment maximum alignment score.
   * @param properRate rate of properly mated pairs.
   * @param properRandomRate rate that mated pairs will appear at random.
   * @param discordantRate rate for discordant mated pairs.
   * @param unmatedRate rate for unmated reads.
   */
  public ReadGroupStats(String id, double meanLength, double fragmentMean, double fragmentStdDev, double gapMean, double gapStdDev, int maxAlignment, double properRate, double properRandomRate, double discordantRate, double unmatedRate) {
    // For testing - set precalculated values
    mLengthDist = null;
    mFragmentSizeDist = null;
    mGapSizeDist = null;
    mMeanLength = meanLength;
    mReadGroup = id;
    mFragmentMean = fragmentMean;
    mFragmentStdDev = fragmentStdDev;
    mGapMean = gapMean;
    mGapStdDev = gapStdDev;
    mMaxAlignment = maxAlignment;
    mProperRate = properRate;
    mProperRandomRate = properRandomRate;
    mDiscordantRate = discordantRate;
    mUnmatedRate = unmatedRate;
    mAlignmentStartIgnored = mMaxAlignment;
  }

  void add(ReadGroupStats stats) {
    mLengthDist.add(stats.mLengthDist);
    mFragmentSizeDist.add(stats.mFragmentSizeDist);
    mGapSizeDist.add(stats.mGapSizeDist);
    mMaxAlignment = Math.max(mMaxAlignment, stats.mMaxAlignment);
    mProper += stats.mProper;
    mDiscordant += stats.mDiscordant;
    mUnmated += stats.mUnmated;
    if (detectOverflow()) {
      throw new NoTalkbackSlimException("Overflow detected in read group statistics calculation for read group: " + id());
    }
  }

  void addLength(int readlen) {
    mLengthDist.add(readlen);
  }
  void addFragmentSize(int value) {
    mFragmentSizeDist.add(value);
  }
  void addGapSize(int value) {
    mGapSizeDist.add(value);
  }
  void addScore(int score) {
    mMaxAlignment = Math.max(mMaxAlignment, score);
  }
  void addUnmated() {
    ++mUnmated;
  }
  void addProper() {
    ++mProper;
  }
  void addDiscordant() {
    ++mDiscordant;
  }
  long getMatedCount() {
    return mProper;
  }
  long getTotalReference() {
    return mTotalReference;
  }
  static String countsHeader() {
    return "#RG-ID"
        + "\ttotal-ref"
        + "\tlength-n"
        + "\tlength-x"
        + "\tlength-x2"
        + "\ttlen-n"
        + "\ttlen-x"
        + "\ttlen-x2"
        + "\tgap-n"
        + "\tgap-x"
        + "\tgap-x2"
        + "\tmax-score"
        + "\tproper"
        + "\tdiscordant"
        + "\tunmated";
  }
  String countsString() {
    return mReadGroup
        + "\t" + mTotalReference
        + "\t" + mLengthDist.count()
        + "\t" + Utils.realFormat(mLengthDist.sum())
        + "\t" + Utils.realFormat(mLengthDist.sumSq())
        + "\t" + mFragmentSizeDist.count()
        + "\t" + Utils.realFormat(mFragmentSizeDist.sum())
        + "\t" + Utils.realFormat(mFragmentSizeDist.sumSq())
        + "\t" + mGapSizeDist.count()
        + "\t" + Utils.realFormat(mGapSizeDist.sum())
        + "\t" + Utils.realFormat(mGapSizeDist.sumSq())
        + "\t" + mMaxAlignment
        + "\t" + mProper
        + "\t" + mDiscordant
        + "\t" + mUnmated;
  }

  void calculate() {
    if (mLengthDist != null) {
      mMeanLength = mLengthDist.count() > 1 ? mLengthDist.mean() : 0;
      mFragmentStdDev = mFragmentSizeDist.count() > 1 ? mFragmentSizeDist.stdDev() + TINY : 0;
      mFragmentMean = mFragmentSizeDist.count() > 1 ? mFragmentSizeDist.mean() : 0;
      mGapStdDev = mGapSizeDist.count() > 1 ? mGapSizeDist.stdDev() + TINY : 0;
      mGapMean = mGapSizeDist.count() > 1 ? mGapSizeDist.mean() : 0;
      mProperRate = (double) (mProper + 1) / mTotalReference / 2.0 + TINY;
      mDiscordantRate = (double) (mDiscordant + 1) / mTotalReference / 2.0 + TINY;
      mUnmatedRate = (double) (mUnmated + 1) / mTotalReference / 2.0 + TINY;

      mProperRandomRate = mDiscordantRate * mUnmatedRate;
      // TODO we think this should be better:
      //mProperRandomRate = mDiscordantRate * (2 * mGapStdDev / mTotalReference);
      mAlignmentStartIgnored = (int) (mMeanLength * ALIGNMENT_START_IGNORED_FRACTION / 100);
    }
  }

  /**
   * A default estimate of the lower bound for a window (inclusive).
   * @return the lower bound.
   */
  public int lo() {
    return -(int) (fragmentMean() + 3 * fragmentStdDev() + meanLength());
  }

  /**
   * A default estimate of the upper bound for a window (exclusive).
   * @return the upper bound.
   */
  public int hi() {
    return (int) (fragmentMean() + 3 * fragmentStdDev()) + 1;
  }

  /**
   * Gets the identifier of the read group
   *
   * @return the read group id
   */
  public String id() {
    return mReadGroup;
  }

  /**
   * Gets the mean of lengths of reads.
   *
   * @return a <code>double</code> value
   */
  public double meanLength() {
    return mMeanLength;
  }

  /**
   * Gets the mean of the absolute mated fragment sizes.
   *
   * @return a <code>double</code> value
   */
  public double fragmentMean() {
    return mFragmentMean;
  }

  /**
   * Gets the standard deviation of the absolute mated fragment sizes.
   *
   * @return a <code>double</code> value
   */
  public double fragmentStdDev() {
    return mFragmentStdDev;
  }

  /**
   * Gets the mean of the absolute mated gap sizes.
   *
   * @return a <code>double</code> value
   */
  public double gapMean() {
    return mGapMean;
  }

  /**
   * Gets the standard deviation of the absolute mated gap sizes.
   *
   * @return a <code>double</code> value
   */
  public double gapStdDev() {
    return mGapStdDev;
  }

  /**
   * Gets the number of bases that may have been incorrectly aligned across a breakpoint due to alignment penalties.
   * @return am <code>int</code> value
   */
  public int alignmentStartIgnored() {
    return mAlignmentStartIgnored;
  }

  /**
   * Gets the maximum alignment score.
   *
   * @return am <code>int</code> value
   */
  public int maxAlignment() {
    return mMaxAlignment;
  }


  /**
   * Gets the rate of proper matings (mate mapping inside
   * expected range) per template base.
   *
   * @return a <code>double</code> value
   */
  public double properRate() {
    return mProperRate;
  }

  /**
   * Gets the rate of proper matings that occur at random
   * (for example, in a deletion region)
   *
   * @return a <code>double</code> value
   */
  public double properRandomRate() {
    return mProperRandomRate;
  }

  /**
   * Gets the rate of discordant matings (mate mapping outside of
   * expected range) per template base.
   *
   * @return a <code>double</code> value
   */
  public double discordantRate() {
    return mDiscordantRate;
  }

  /**
   * Gets the rate of unmated reads per template base.
   *
   * @return a <code>double</code> value
   */
  public double unmatedRate() {
    return mUnmatedRate;
  }

  @Override
  public String toString() {
    return id()
        + "\t" + Utils.realFormat(meanLength(), 5)
        + "\t" + Utils.realFormat(fragmentMean(), 5)
        + "\t" + Utils.realFormat(fragmentStdDev(), 5)
        + "\t" + Utils.realFormat(gapMean(), 5)
        + "\t" + Utils.realFormat(gapStdDev(), 5)
        + "\t" + Utils.realFormat(maxAlignment(), 5)
        + "\t" + Utils.realFormat(properRate(), 6)
        + "\t" + Utils.realFormat(properRandomRate(), 6)
        + "\t" + Utils.realFormat(discordantRate(), 6)
        + "\t" + Utils.realFormat(unmatedRate(), 6);
  }

  boolean detectOverflow() {
    return mProper < 0 || mDiscordant < 0 || mUnmatedRate < 0 || mLengthDist.count() < 0;
  }

  /**
   * Check if this statistics object contains valid values.
   * @return <code>true</code> if valid, <code>false</code> otherwise
   */
  public boolean isValid() {
    return checkDoubleNonNegativeNatural(mMeanLength)
        && checkDoubleNonNegativeNatural(mFragmentMean)
        && checkDoubleNonNegativeNatural(mFragmentStdDev)
        && checkDoubleNatural(mGapMean)
        && checkDoubleNonNegativeNatural(mGapStdDev)
        && checkDoubleNonNegativeNatural(mProperRate)
        && checkDoubleNonNegativeNatural(mProperRandomRate)
        && checkDoubleNonNegativeNatural(mDiscordantRate)
        && checkDoubleNonNegativeNatural(mUnmatedRate);
  }

  private boolean checkDoubleNonNegativeNatural(double value) {
    return checkDoubleNatural(value) && value >= 0.0;
  }

  private boolean checkDoubleNatural(double value) {
    return !Double.isInfinite(value) && !Double.isNaN(value);
  }

  /**
   * Test code from command line: print out stats for supplied read group stats files.
   *
   * @param args a <code>String</code> value
   * @exception Exception if an error occurs.
   */
  public static void main(String[] args) throws Exception {
    final File[] files = new File[args.length];
    for (int i = 0; i < args.length; ++i) {
      files[i] = new File(args[i]);
    }
    loadReadGroupStats(null, files);
  }
}
