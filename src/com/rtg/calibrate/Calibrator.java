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

package com.rtg.calibrate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SuperCigarValidator;
import com.rtg.util.Environment;
import com.rtg.util.Histogram;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.AsynchInputStream;

import htsjdk.samtools.SAMRecord;

/**
 * Measures the actual machine error rates observed in SAM files.
 *
 * Maintains a Hypercube of statistics objects,
 * indexed by a vector of Covariate objects.
 *
 */
public class Calibrator {

  private static final String CALIBRATE_VERSION = "calibrate v3.0";
  protected static final String HEADER_CHAR = "@";
  protected static final String COVAR = HEADER_CHAR + "covar";
  /** Constant for retrieving distribution */
  public static final String CGOVER_DIST = HEADER_CHAR + "cgover";
  /** Constant for retrieving distribution */
  public static final String CGGAP_DIST = HEADER_CHAR + "cggap";
  /** Constant for retrieving distribution */
  public static final String DEL_DIST = HEADER_CHAR + "del";
  /** Constant for retrieving distribution */
  public static final String INS_DIST = HEADER_CHAR + "ins";
  /** Constant for retrieving distribution */
  public static final String MNP_DIST = HEADER_CHAR + "mnp";
  /** Constant for retrieving distribution */
  public static final String NH_DIST = HEADER_CHAR + "nh";

  private static final String STATS_COLUMNS = "\tequal\tdiff\tins\tdel";
  private static final int NUM_STATS_COLUMNS = 4;
  private static final String REFERENCE_SIZE = HEADER_CHAR + "sequence";

  /** The index dimensions of the statistics hypercube. */
  protected final Covariate[] mCovariates;

  /** The actual statistics collection. */
  protected CalibrationStats[] mStats;
  protected final CalibratorCigarParser mParser;
  private String mTemplateName = null;
  protected byte[] mTemplate;
  protected int mTemplateLength;
  protected byte[] mRead;
  protected BitSet mTemplateMask = null;
  private byte[] mCgQualities = null;
  protected SAMRecord mSamRec;
  protected String mReadGroup;
  final ReferenceRegions mRegions;
  protected Map<String, Histogram> mDistributions = new HashMap<>();
  private final Map<String, Integer> mSequenceLengths = new HashMap<>();

  /**
   * Create a calibrator where statistics will only be accumulated within the specified reference regions
   * @param vars the sequence of covariates that may influence the error rates/qualities.
   * @param regions regions within which calibration statistics should be accumulated.
   */
  public Calibrator(Covariate[] vars, ReferenceRegions regions) {
    mCovariates = vars;
    mStats = new CalibrationStats[maxSize(vars)];
    mParser = new CalibratorCigarParser(this);
    mRegions = regions;
    if (mRegions != null) {
      mSequenceLengths.putAll(mRegions.coveredLengths());
    }
  }

  /**
   * Construct a map from sequence name to length covered by BED regions.
   * @param reader the sequences reader to construct the map from
   * @param regions only measure the parts within these regions.
   * @return the map from sequence name to covered length
   * @throws java.io.IOException if your sequences reader is dodgy.
   */
  public static Map<String, Integer> getSequenceLengthMap(SequencesReader reader, ReferenceRegions regions) throws IOException {
    final Map<String, Integer> lengthMap = new HashMap<>();
    for (int i = 0; i < reader.numberSequences(); ++i) {
      lengthMap.put(reader.name(i), regions.coveredLength(reader.name(i)));
    }
    return lengthMap;
  }

  /**
   * Construct a map from sequence name to non-n length calculated from the reference.
   * @param reader the sequences reader to construct the map from
   * @param restriction if not null, only include information for the sequence contained in the restriction
   * @return the map from sequence name to sequence id
   * @throws java.io.IOException if your sequences reader is dodgy.
   */
  public static Map<String, Integer> getNonNSequenceLengthMap(SequencesReader reader, RegionRestriction restriction) throws IOException {
    final Map<String, Integer> lengthMap = new HashMap<>();
    for (int i = 0; i < reader.numberSequences(); ++i) {
      final String name = reader.name(i);
      if (restriction != null && !name.equals(restriction.getSequenceName())) {
        continue;
      }

      int length = reader.length(i);

      // Adjust lengths to subtract off n's - this ensures that the calibrated coverage is more correct
      // for chromosomes with large chunks of n's, since we won't have mapping coverage
      final byte[] refNts = new byte[length];
      reader.read(i, refNts);
      for (final byte b : refNts) {
        if (b == DnaUtils.UNKNOWN_RESIDUE) {
          --length;
        }
      }
      Diagnostic.developerLog("Length of sequence " + name + " for calibration is " + length + " (" + refNts.length + " raw)");

      lengthMap.put(name, length);
    }
    return lengthMap;
  }

  /**
   * Construct a map from sequence name to non-n length calculated from the reference.
   * @param reader the sequences reader to construct the map from
   * @param restriction if not null, only include information for the sequence contained in the restriction
   * @return the map from sequence name to sequence id
   * @throws java.io.IOException if your sequences reader is dodgy.
   */
  public static ReferenceRegions getNonNRegions(SequencesReader reader, RegionRestriction restriction) throws IOException {
    final ReferenceRegions r = new ReferenceRegions();
    for (int i = 0; i < reader.numberSequences(); ++i) {
      final String name = reader.name(i);
      if (restriction != null && !name.equals(restriction.getSequenceName())) {
        continue;
      }
      final int length = reader.length(i);
      final byte[] refNts = new byte[length];
      reader.read(i, refNts);
      int start = 0;
      int end = 0;
      for (int i1 = 0; i1 < refNts.length; ++i1) {
        if (refNts[i1] != DnaUtils.UNKNOWN_RESIDUE) {
          end = i1 + 1;
        } else {
          if (start < end) {
            r.add(name, start, end);
          }
          start = i1 + 1;
        }
      }
      if (start < end) {
        r.add(name, start, end);
      }
      Diagnostic.developerLog("Length of sequence " + name + " for calibration is " + r.coveredLength(name) + " (" + refNts.length + " raw) ");
    }
    return r;
  }

  /**
   * Initialize calibration.
   * @param calibrationFiles calibration files
   * @return the calibrator
   * @exception java.io.IOException if an I/O error occurs.
   */
  public static Calibrator initCalibrator(Collection<File> calibrationFiles) throws IOException {
    Covariate[] covariates = null;
    if (calibrationFiles.size() > 0) {
      for (final File f : calibrationFiles) {
        covariates = getCovariateSet(f);
        break;
      }
    } else {
      return null;
    }
    if (covariates == null) {
      throw new IllegalStateException("no calibration covariates found");
    }
    final Calibrator c = new Calibrator(covariates, null);
    for (final File f : calibrationFiles) {
      c.accumulate(f);
    }
    return c;
  }

  CalibratorCigarParser getParser() {
    return mParser;
  }

  /**
   * @param vars the index dimensions of the hypercube.
   * @return the maximum number of statistics entries in the hypercube.
   */
  protected static int maxSize(Covariate[] vars) {
    long maxSize = 1;
    final int maxValue = Integer.MAX_VALUE;
    for (final Covariate var : vars) {
      maxSize *= var.newSize();
      if (maxSize > maxValue) {
        throw new SlimException("Too many covariates during calibration");
      }
    }
    return (int) maxSize;
  }

  /**
   * Set the sequence lengths map
   * @param sequenceLengths map of sequence name to length
   */
  public void setSequenceLengths(Map<String, Integer> sequenceLengths) {
    mSequenceLengths.putAll(sequenceLengths);
  }

  public Map<String, Integer> getSequenceLengths() {
    return mSequenceLengths;
  }

  /**
   * The recalibrate tool may be restricted to the regions within a bed file, If so the total lengths of these regions
   * will be stored in this class and this function will return true.
   * @return true if this calibrator contains sequence length information
   */
  public boolean hasLengths() {
    return mSequenceLengths.size() > 0;
  }
  
  boolean inRange(int templatePosition) {
    return mTemplateMask == null || mTemplateMask.get(templatePosition);
  }

  /**
   * Get a histogram.
   * @param label label for the histogram
   * @return the histogram, or null if doesn't exist
   */
  Histogram getHistogram(String label) {
    return mDistributions.get(label);
  }

  /**
   * Get histogram, or create if doesn't exist.
   * Not thread safe - check first with <code>hasHistogram</code>.
   * @param label label for the histogram
   * @return the histogram
   */
  Histogram getOrCreateHistogram(String label) {
    final Histogram ret = mDistributions.get(label);
    if (ret != null) {
      return ret;
    }
    final Histogram ret1 = new Histogram();
    mDistributions.put(label, ret1);
    return ret1;
  }

  /**
   * Gets a histogram. Pretty please don't modify it.
   * @param type distribution as defined above
   * @param readGroup read group id
   * @return histogram containing distribution counts
   */
  public Histogram getHistogram(String type, String readGroup) {
    return getHistogram(type + ":" + readGroup);
  }

  /**
   * Find if a given histogram already exists.
   *
   * @param dist a histogram type
   * @param readGroup name of read group
   * @return true if it already has the desired histogram.
   */
  public boolean hasHistogram(String dist, String readGroup) {
    return mDistributions.containsKey(dist + ":" + readGroup);
  }

  /**
   * @param index the index of the covariate as given by {@link Calibrator#getCovariateIndex(com.rtg.calibrate.CovariateEnum)}
   * @return the covariate object
   */
  public Covariate getCovariate(int index) {
    return mCovariates[index];
  }

  /**
   * Get a blank covariate set from the given calibration file
   * @param cal file containing calibration data
   * @return the covariates
   * @throws IOException if an IO error occurs
   */
  public static Covariate[] getCovariateSet(File cal) throws IOException {
    try (BufferedReader input = new BufferedReader(new FileReader(cal))) {
      String line;
      while ((line = input.readLine()) != null) {
        if (line.startsWith("#")) {
          continue;
        }
        if (line.startsWith(COVAR)) {
          final String[] names = line.split("\t");
          if (!line.endsWith(STATS_COLUMNS)) {
            throw new NoTalkbackSlimException("calibration file \"" + cal.getPath() + "\" does not have correct stats columns: " + line);
          }
          final int start = 1;
          final int end = names.length - NUM_STATS_COLUMNS;
          final Covariate[] ret = new Covariate[end - start];
          for (int i = start; i < end; ++i) {
            final String[] cvsplit = names[i].split(":");
            final String cvname = cvsplit[0];
            final int length = cvsplit.length > 1 ? Integer.parseInt(cvsplit[1]) : 0;
            ret[i - start] = CovariateEnum.getCovariate(null, CovariateEnum.valueOf(cvname.toUpperCase(Locale.ROOT)), length);
          }
          return ret;
        }
      }
    }
    throw new IOException("File: " + cal + " did not contain covariate definitions");
  }

  /**
   * Reads another calibration file and adds its statistics into this one.
   * Not thread safe (due to <code>getOrCreateHistogram</code> calls)
   *
   * @param file path to read from.
   * @throws IOException if read or input format errors occur.
   */
  public void accumulate(File file) throws IOException {
    try (InputStream is = new FileInputStream(file)) {
      accumulate(is, file.getPath());
    }
  }

  void accumulate(InputStream inputStream, String path) throws IOException {
    try (BufferedReader input = new BufferedReader(new InputStreamReader(new AsynchInputStream(inputStream)))) {
      String line;
      while ((line = input.readLine()) != null) {
        if (line.startsWith("#")) {
          continue;
        }
        final String[] names = StringUtils.split(line, '\t');
        if (names[0].startsWith(REFERENCE_SIZE)) {
          final int secondTab = names[0].length() + names[1].length() + 2;
          try {
            mSequenceLengths.put(line.substring(secondTab), Integer.parseInt(names[1]));
          } catch (final NumberFormatException e) {
            throw new NoTalkbackSlimException("calibration file had an invalid sequence length: " + names[1]);
          }

        } else if (names[0].startsWith(NH_DIST) || names[0].startsWith(MNP_DIST) || names[0].startsWith(INS_DIST)
            || names[0].startsWith(DEL_DIST) || names[0].startsWith(CGGAP_DIST) || names[0].startsWith(CGOVER_DIST)) {
          final Histogram hist = getOrCreateHistogram(names[0]);
          if (names.length > 1) {
            hist.addHistogram(line.substring(names[0].length() + 1));
          }
        } else if (names[0].startsWith(COVAR)) {
          if (!line.endsWith(STATS_COLUMNS)) {
            throw new NoTalkbackSlimException("calibration file \"" + path + "\" does not have correct stats columns: " + line);
          }
          if ((names.length - NUM_STATS_COLUMNS - 1) != mCovariates.length) {
            throw new NoTalkbackSlimException("calibration file \"" + path + "\" has mismatching covariates: " + line);
          }
          for (int i = 0; i < mCovariates.length; ++i) {
            final String[] nameParts = StringUtils.split(names[i + 1], ':');
            final String[] covNameParts = StringUtils.split(mCovariates[i].name(), ':');
            if (!nameParts[0].equals(covNameParts[0])) {
              throw new NoTalkbackSlimException("calibration file \"" + path + "\" contains unexpected covariate: " + nameParts[0]);
            }
          }
          break;
        } else {
          throw new NoTalkbackSlimException("calibration file \"" + path + "\" has bad header line: " + line);
        }
      }

      while ((line = input.readLine()) != null) {
        final String[] field = StringUtils.split(line, '\t');
        if (field.length != mCovariates.length + NUM_STATS_COLUMNS) {
          throw new NoTalkbackSlimException("calibration file \"" + path + "\" contains invalid line: " + line);
        }
        int pos = 0;
        int i = 0;
        boolean mustResize = false;
        final int[] values = new int[mCovariates.length];
        for (; i < mCovariates.length; ++i) {
          final int val = mCovariates[i].parse(field[i]);
          values[i] = val;
          mustResize |= mCovariates[i].sizeChanged();
          pos = pos * mCovariates[i].newSize() + val;
        }
        if (mustResize) {
          expandStats();
        }
        if (mStats[pos] != null) {
          assert Arrays.equals(values, mStats[pos].getValues());
        } else {
          mStats[pos] = new CalibrationStats(values);
        }
        final CalibrationStats stats = mStats[pos];
        stats.seenEquals(Long.parseLong(field[i++]));
        stats.seenMnp(Long.parseLong(field[i++]));
        stats.seenInsert(Long.parseLong(field[i++]));
        stats.seenDelete(Long.parseLong(field[i++]));
        assert i == field.length;
      }
    }
  }

  /**
   * Accumulate the results of another compatible calibrator into this calibrator.
   * @param cal other calibrator
   */
  public void accumulate(Calibrator cal) {
    // Check covariates match
    if (mCovariates.length != cal.mCovariates.length) {
      throw new RuntimeException("Missing covariates");
    }
    for (int k = 0; k < mCovariates.length; ++k) {
      if (!getCovariate(k).name().equals(cal.getCovariate(k).name())) {
        throw new RuntimeException("Covariates mismatch");
      }
    }
    // Merge hypercubes
    assert mStats.length == cal.mStats.length;
    for (int k = 0; k < mStats.length; ++k) {
      if (mStats[k] != null) {
        if (cal.mStats[k] != null) {
          mStats[k].accumulate(cal.mStats[k]);
        }
      } else {
        mStats[k] = cal.mStats[k];
      }
    }
    // Merge histograms
    for (final Map.Entry<String, Histogram> e : mDistributions.entrySet()) {
      final Histogram h = cal.getHistogram(e.getKey());
      if (h != null) {
        e.getValue().addHistogram(h);
      }
    }
    for (final Map.Entry<String, Histogram> e : cal.mDistributions.entrySet()) {
      final Histogram h = getHistogram(e.getKey());
      if (h == null) {
        mDistributions.put(e.getKey(), e.getValue());
      }
    }
  }

  /**
   * @param ce type of covariate
   * @return the index of the covariate for use in other methods,
   *         or -1 if the requested covariate is not in the hypercube.
   */
  public int getCovariateIndex(CovariateEnum ce) {
    for (int i = 0; i < mCovariates.length; ++i) {
      if (mCovariates[i].getType() == ce) {
        return i;
      }
    }
    return -1;
  }

  /** Expands the statistics array and recalculates all the indices. */
  private void expandStats() {
    final int start = mStats.length - 1;
    mStats = Arrays.copyOf(mStats, maxSize(mCovariates));
    for (int i = start; i > 0; --i) {
      if (mStats[i] != null) {
        final int newPos = mStats[i].getNewIndex(mCovariates);
        if (newPos != i) {
          mStats[newPos] = mStats[i];
          mStats[i] = null;
        }
      }
    }
    for (final Covariate cov : mCovariates) {
      cov.resized();
    }
  }

  protected void writeHistogram(String name, Histogram hist, BufferedWriter out) throws IOException {
    if (hist.getLength() > 0) {
      out.write(name + "\t" + hist.toString() + StringUtils.LS);
    } else {
      out.write(name + StringUtils.LS);
    }
  }

  /**
   * Get the sum of all statistics for a given covariate value. Or all statistics
   * if given covariate is not in use.
   * @param covariate covariate type. null means all
   * @param covariateValue value of covariate
   * @return sum over all bins with given covariate value
   */
  public CalibrationStats getSums(CovariateEnum covariate, String covariateValue) {
    int covNum = -1;
    int covVal = -1;
    for (int i = 0; i < mCovariates.length; ++i) {
      if (mCovariates[i].getType() == covariate) {
        covVal = mCovariates[i].parse(covariateValue);
        covNum = i;
        break;
      }
    }

    final CalibrationStats newStats = new CalibrationStats(new int[] {covVal});
    for (final CalibrationStats mStat : mStats) {
      if (mStat != null) {
        if (covNum == -1 || mStat.getCovariateValue(covNum) == covVal) {
          newStats.accumulate(mStat);
        }
      }
    }
    return newStats;
  }

  /**
   * Iterates through the subset of the hypercube that satisfies <code>query</code>
   * and calls the call-back method <code>proc.process</code> for each one.
   *
   * @param proc the call-back method to call on every statistics object that satisfies the query.
   * @param query specifies which dimensions have a fixed value and which to range over.
   */
  public void processStats(StatsProcessor proc, QuerySpec query) {
    processStats(proc, query, new int[mCovariates.length], 0, 0);
  }

  /** recursive helper method for <code>processStats</code>. */
  private void processStats(StatsProcessor proc, QuerySpec query, int[] values, int position, int covariateIndex) {
    if (covariateIndex == mCovariates.length) {
      proc.process(values, mStats[position]);
      return;
    }
    final int size = mCovariates[covariateIndex].size();
    if (query.mCovariateValues[covariateIndex] >= 0) {
      values[covariateIndex] = query.mCovariateValues[covariateIndex];
      processStats(proc, query, values, position * size + query.mCovariateValues[covariateIndex], covariateIndex + 1);
    } else {
      for (int i = 0; i < size; ++i) {
        values[covariateIndex] = i;
        processStats(proc, query, values, position * size + i, covariateIndex + 1);
      }
    }
  }

  /**
   * Dumps all the statistics out into a text file.
   *
   * @param file The file to write to.
   * @throws IOException if an IO error occurs during output.
   */
  public void writeToFile(File file) throws IOException {
    try (BufferedWriter out = new BufferedWriter(new FileWriter(file))) {
      write(out);
    }
  }

  /**
   * Dumps all the statistics out into an output stream.
   * @param stream the output stream to write to.
   * @throws IOException if an IO error occurs during output.
   */
  public void writeToStream(OutputStream stream) throws IOException {
    try (BufferedWriter out = new BufferedWriter(new OutputStreamWriter(stream))) {
      write(out);
    }
  }

  /**
   * Dumps all the statistics out into a buffered writer.
   * @param out the buffered writer
   * @throws IOException if an IO error occurs during output.
   */
  private void write(BufferedWriter out) throws IOException {
    out.write("#Version " + Environment.getVersion() + ", " + CALIBRATE_VERSION + StringUtils.LS);
    out.write("#CL\t" + CommandLine.getCommandLine() + StringUtils.LS);
    final Set<String> histogramLabels = new TreeSet<>(mDistributions.keySet());
    for (final String label : histogramLabels) {
      writeHistogram(label, getHistogram(label), out);
    }
    for (final String sequence : new TreeSet<>(mSequenceLengths.keySet())) {
      out.write(REFERENCE_SIZE + "\t" + mSequenceLengths.get(sequence) + "\t" + sequence + StringUtils.LS);
    }

    out.write(COVAR + "\t" + toString() + STATS_COLUMNS + StringUtils.LS);  // tab separated
    for (final CalibrationStats stats : mStats) {
      if (stats != null) {
        out.write(stats.outputString(mCovariates) + StringUtils.LS);
      }
    }
  }

  /**
   * create a new blank query
   * @return the query
   */
  public QuerySpec initQuery() {
    return new QuerySpec();
  }

  protected CalibrationStats findStats(CalibratorCigarParser currPos) throws BadSuperCigarException {
    int pos = 0;
    boolean mustResize = false;
    final int[] values = new int[mCovariates.length];
    for (int i = 0; i < mCovariates.length; ++i) {
      final Covariate var = mCovariates[i];
      final int val = var.value(mSamRec, currPos);
      values[i] = val;
      mustResize |= var.sizeChanged();
      pos = pos * var.newSize() + val;
    }
    if (mustResize) {
      expandStats();
    }
    CalibrationStats stats = mStats[pos];
    if (stats == null) {
      stats = new CalibrationStats(values);
      mStats[pos] = stats;
    }
    return stats;
  }

  /**
   * Must call this each time we start a new template sequence.
   * @param name the reference sequence name
   * @param template template nucleotides
   * @param length amount of array used by template
   */
  public void setTemplate(String name, byte[] template, int length) {
    setTemplate(name, 0, template, length);
  }

  /**
   * Must call this each time we start a new template sequence. This method allows
   * supplying a template buffer that covers only a portion of the template sequence.
   * @param name the reference sequence name
   * @param templateOffset the zero-based coordinate that the first byte of the provided template corresponds to
   * @param template template nucleotides
   * @param length amount of array used by template
   */
  public void setTemplate(String name, int templateOffset, byte[] template, int length) {
    if (mRegions != null && !name.equals(mTemplateName)) {
      mTemplateMask = mRegions.mask(name);
    }
    mTemplateName = name;
    mTemplate = template;
    mTemplateLength = length;
    getParser().setTemplate(templateOffset, mTemplate, mTemplateLength);
  }

  /**
   * @return a tab-separated sequence of the covariate names.
   */
  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    String sep = "";
    for (final Covariate var : mCovariates) {
      sb.append(sep).append(var.name());
      sep = StringUtils.TAB;
    }
    return sb.toString();
  }

  private Histogram mNhHistogram;

  /**
   * Records error statistics for all positions along the given read.
   * It parses the cigar (or super cigar) and accumulates counts of all events.
   * Not thread safe (due to <code>getOrCreateHistogram</code> calls)
   *
   * @param sam The SAM record for a read.
   */
  public void processRead(SAMRecord sam) {
    mSamRec = sam;
    final String readGroup = ReadGroupUtils.getReadGroup(sam);
    if (!readGroup.equals(mReadGroup)) {
      mReadGroup = readGroup;
      mNhHistogram = getOrCreateHistogram(NH_DIST + ":" + mReadGroup);
    }
    if (ReadGroupUtils.UNKNOWN_RG.equals(mReadGroup)) {
      throw new NoTalkbackSlimException("quality calibration requires a read group to be specified");
    }
    final int zeroBasedStart = sam.getAlignmentStart() - 1;
    final int zeroBasedEnd = sam.getAlignmentEnd();
    if (mRegions != null && !mRegions.overlapped(sam.getReferenceName(), zeroBasedStart, zeroBasedEnd)) {
      return;
    }
    final Integer inh = SamUtils.getNHOrIH(sam);
    final int nh = inh == null ? 1 : inh;
    mNhHistogram.increment(nh);
    if (nh > 1) {
      return;
    }
    final int mapq = sam.getMappingQuality();
    if (mapq == 0) {
      return;
    }
    getParser().setTemplateStart(zeroBasedStart);
    try {
      final byte[] baseQualities = sam.getBaseQualities();
      if (sam.hasAttribute(SamUtils.CG_SUPER_CIGAR)) {
        final String superCigar = sam.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
        final String readDelta = sam.getStringAttribute(SamUtils.CG_READ_DELTA);
        getParser().setCigar(superCigar, readDelta);

        final int flags = sam.getFlags();
        final boolean first = (flags & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
        final boolean rc = (flags & SamBamConstants.SAM_READ_IS_REVERSE) != 0;
        final String xqField = sam.getStringAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY);
        if (baseQualities != null && baseQualities.length > 0) {
          final int xqlength = xqField == null ? 0 : xqField.length();
          if (mCgQualities == null || baseQualities.length + xqlength != mCgQualities.length) {
            mCgQualities = new byte[baseQualities.length + xqlength];
          }
          SuperCigarValidator.unrollSuperCigarQualities(mCgQualities, baseQualities, xqField, first, rc);
          getParser().setQualities(mCgQualities);
        } else {
          getParser().setQualities(null);
        }
      } else {
        final String cigar = sam.getCigarString();
        final int readLength = sam.getReadLength();
        if (mRead == null || readLength > mRead.length) {
          mRead = new byte[readLength];
        }
        System.arraycopy(sam.getReadBases(), 0, mRead, 0, readLength);
        mRead = DnaUtils.encodeArray(mRead, readLength);
        getParser().setStandardCigar(cigar, mRead, readLength);
        if (baseQualities != null && baseQualities.length > 0) {
          getParser().setQualities(baseQualities);
        } else {
          getParser().setQualities(null);
        }
      }
      getParser().parse();
    } catch (final BadSuperCigarException e) {
      Diagnostic.developerLog("Ignored SAM record due to " + e.getMessage() + " query=" + sam.getReadName()
          + " ref=" + sam.getReferenceName()
          + ":" + sam.getAlignmentStart());
    }
    mSamRec = null;
  }

  /**
   * Specifies fixed values for some dimensions of the hypercube and allows the other dimensions to vary.
   */
  public class QuerySpec {
    private final int[] mCovariateValues = new int[mCovariates.length];

    /**
     * Constructor.  Allows all dimensions to vary.
     */
    public QuerySpec() {
      Arrays.fill(mCovariateValues, -1);
    }

    /**
     * Restrict the given covariate to have the given value.
     * @param cov which covariate to restrict.
     * @param covariateValue the value that covariate must have.
     * @return true iff <code>cov</code> is one of the dimensions of the current hypercube.
     */
    public boolean setValue(CovariateEnum cov, int covariateValue) {
      final int i = getCovariateIndex(cov);
      if (i != -1) {
        mCovariateValues[i] = covariateValue;
        return true;
      }
      return false;
    }
  }
}
