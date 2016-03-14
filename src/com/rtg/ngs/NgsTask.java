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
package com.rtg.ngs;


import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.ChrStats;
import com.rtg.calibrate.Recalibrate;
import com.rtg.index.Index;
import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.NgsHashLoop;
import com.rtg.index.hash.ngs.NgsHashLoopImpl;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.index.hash.ngs.ReadCallImplementation;
import com.rtg.index.hash.ngs.ReadEncoder;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.TemplateCallImplementation;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.ParamsTask;
import com.rtg.ngs.longread.LongReadTask;
import com.rtg.position.output.PositionParams;
import com.rtg.reader.CgUtils;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderUtils;
import com.rtg.report.MapFSummaryReport;
import com.rtg.report.MapReport;
import com.rtg.report.MapSummaryReport;
import com.rtg.report.ReportType;
import com.rtg.usage.UsageMetric;
import com.rtg.util.Environment;
import com.rtg.util.MathUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;

/**
 * Takes reads and a template and generates mappings.
 */
public class NgsTask extends ParamsTask<NgsParams, MapStatistics> {

  /**
   * Construct a new build and search.
   * @param params parameters for the build and search.
   * @param defaultOutput place to send output.
   * @param metric container for number of reads processed for usage tracking purposes
   */
  public NgsTask(final NgsParams params, final OutputStream defaultOutput, final UsageMetric metric) {
    super(params, defaultOutput, createStatistics(params), metric);
    if (params.paired()) {
      if (params.buildFirstParams().numberSequences() != params.buildSecondParams().numberSequences()) {
        throw new RuntimeException("Number of reads in first and second read sets must be equal.");
      }
      if (!params.useLongReadMapping()) {
        if (params.buildFirstParams().maxLength() != params.buildSecondParams().maxLength()) {
          throw new RuntimeException("Length of reads in first and second read sets must be equal.");
        }
      }
    }
  }

  private static MapStatistics createStatistics(final NgsParams params) {
    if (params.outputParams().filter().outputFilter() == OutputFilter.SAM_UNFILTERED) {
      if (params.paired()) {
        return new MapFilterPairedMapStatistics(params.directory());
      } else {
        return new MapFilterSingleEndMapStatistics(params.directory());
      }
    } else if (params.paired()) {
      return new PairedEndMapStatistics(params.outputParams().outputUnmated(), params.directory());
    } else {
      return new SingleEndMapStatistics(params.directory());
    }
  }

  /**
   * @return the parameters for this task.
   */
  @Override
  public NgsParams parameters() {
    return mParams;
  }

  /** runs the task */
  @Override
  protected void exec() throws IOException {
    Diagnostic.developerLog("NGSParams" + LS + mParams.toString());
    //make all the components we need
    final OneShotTimer fullTimer = new OneShotTimer("total_time");
    assert mParams.searchParams().numberSequences() < Integer.MAX_VALUE : mParams.buildFirstParams().numberSequences();
    if (mParams.useLongReadMapping()) {
      buildQueryLongRead(mParams, mStatistics, mUsageMetric);
    } else {
      indexThenSearchShortReads(mParams, mStatistics, mUsageMetric);
    }
    logMemStats("Free memory pre-GC ");
    System.gc();
    logMemStats("Free memory post-GC ");
    fullTimer.stopLog();
    chrStatsCheck();
    mapReport();
  }

  /** Put free-memory statistics into the log, if we can. */
  private static void logMemStats(final String when) {
    try {
      Diagnostic.developerLog(when + Environment.getFreeMemory());
    } catch (final IllegalStateException e) {
      Diagnostic.developerLog("Cannot get free memory statistics");
    }
  }

  private List<File> findCalibrationFiles(final File dir) throws IOException {
    return Arrays.asList(FileUtils.listFiles(dir, new FilenameFilter() {
      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(Recalibrate.EXTENSION);
      }
    }));
  }

  private void chrStatsCheck() throws IOException {
    if (mParams.outputParams().calibrate() && mParams.sex() != null) {
      // This checking assumes presence of calibration
      final Collection<File> calibrationFiles = findCalibrationFiles(mParams.outputParams().directory());
      if (!calibrationFiles.isEmpty()) {
        final Calibrator c = Calibrator.initCalibrator(calibrationFiles);
        if (c != null) {
          final String sample = mParams.outputParams().readGroup().getSample();
          final Map<String, String> readGroupToSampleId = Collections.singletonMap(mParams.outputParams().readGroup().getReadGroupId(), sample);
          final Map<String, Integer> sequenceLengthMap = c.hasLengths() ? c.getSequenceLengths() : Calibrator.getNonNSequenceLengthMap(mParams.searchParams().reader(), (RegionRestriction) null);
          final CalibratedPerSequenceExpectedCoverage expectedCoverages = new CalibratedPerSequenceExpectedCoverage(c, sequenceLengthMap, readGroupToSampleId, null);
          final ChrStats cc = new ChrStats(mParams.searchParams().reader());
          if (cc.referenceOk()) {
            final ChrStats.ChrStatsResult res = cc.runCheckAndReport(expectedCoverages, sample, mParams.sex());
            if (res.getObservedSex() != null && res.getClaimedSex() != res.getObservedSex()) {
              Diagnostic.warning(sample + " specified " + ChrStats.sexString(res.getClaimedSex()) + " but appears to be " + ChrStats.sexString(res.getObservedSex()));
            }
          }
        }
      }
    }
  }

  /**
   *   TODO This is fairly ugly (specific to map report, but tidying up the API to allow this to be more module
   *   specific is currently more effort than is worthwhile)
   */
  private void mapReport() throws IOException {
    final MapSummaryReport mr;
    if (new File(mParams.outputParams().directory(), MapReportData.MAP_REPORT_FILE_NAME).exists()) {
      mr = new MapReport(null); //this filter params ought to be ignored, since the data exists...
    } else if (mParams.outputParams().sdf()) {
      // Probably mapf
      mr = new MapFSummaryReport();
      mr.setTitle("Filter Report");

    } else {
        mr = new MapSummaryReport();
    }
    mr.setCommandLines(Collections.singletonList(CommandLine.getCommandLine()));
    mr.setParams(mParams);
    final List<File> samples = new ArrayList<>();
    final File left = mParams.buildFirstParams().reader().path();
    final File right;
    if (mParams.buildSecondParams() != null) {
      right = mParams.buildSecondParams().reader().path();
    } else {
      right = null;
    }
    if (right != null) {
      if (ReaderUtils.isSDF(left) && ReaderUtils.isSDF(right) && left.getParent().equals(right.getParent())
          && left.getName().equals("left")
          && right.getName().equals("right")) {
        samples.add(new File(left.getParent()));
      } else {
        samples.add(left);
        samples.add(right);
      }
    } else {
      samples.add(left);
    }
    mr.setSampleFiles(samples);
    mr.generateReport(ReportType.HTML, mParams.directory(), mParams.directory());
  }

  /**
   * Build and search long read code
   * @param params {@link NgsParams} object
   * @param statistic mapping statistics at the end
   * @param usageMetric accumulates the total number of nucleotides read.
   * @throws IOException if error
   */
  public static void buildQueryLongRead(final NgsParams params, final MapStatistics statistic, final UsageMetric usageMetric) throws IOException {
    final PositionParams posParams = params.toPositionParams();
    if (!checkReadCount(params)) {
      throw new SlimException("Read dataset too large, try running in multiple smaller chunks using --start-read and --end-read parameters");
    }

    final Index index = LongReadTask.build(posParams, usageMetric, params.indexFilter());
    final OutputFilter filter = params.outputParams().outFilter();
    try (OutputProcessor outProcessor = filter.makeProcessor(params, statistic)) {
      LongReadTask.search(posParams, outProcessor, index);
      outProcessor.finish();
    }
  }

  /**
   * Returns true if number of read arms is less than {@link Integer#MAX_VALUE}
   * @param params the params to check
   * @return true if number is within bounds, false otherwise
   */
  private static boolean checkReadCount(NgsParams params) {
    long numReads = params.buildFirstParams().numberSequences();
    if (params.paired()) {
      numReads <<= 1;
    }
    return numReads <= Integer.MAX_VALUE;
  }

  private static void indexThenSearchShortReads(final NgsParams params, final MapStatistics statistics, final UsageMetric usageMetric) throws IOException {
    final long pMask = 0x1FFFFL;
    final Integer numberThreads = params.numberThreads();
    final int threadBits = MathUtils.ceilPowerOf2Bits(numberThreads - 1);
    final HashFunctionFactory factory = params.maskParams().maskFactory((int) params.getMaxReadLength());
    final long numSeqs = params.buildFirstParams().numberSequences() + (params.paired() ? params.buildSecondParams().numberSequences() : 0);
    final CreateParams indexParams = new CreateParams(numSeqs, factory.hashBits(), factory.windowBits(), NgsParams.calculateValueBitsShortReads(params.buildFirstParams().numberSequences(), params.paired()), params.compressHashes(), true, false, false);
    Diagnostic.developerLog("Index params: " + indexParams.toString());
    final NgsHashLoopImpl hashLoop = new NgsHashLoopImpl(params.buildFirstParams().numberSequences(), params.outputParams().progress(), 0x3FFFFL, ((pMask + 1L) << threadBits) - 1L);
    hashLoop.setThreadPadding(params.calculateThreadPadding());
    usageMetric.setMetric(indexThenSearchShortReads(params, hashLoop, statistics, indexParams));
  }

  private static final int INDEX_USAGE_REPORTING_THRESHOLD = 50;

  /**
   * Run query on short read data
   * @param params {@link NgsParams} object
   * @param shl hash loop
   * @param statistics mapping statistics
   * @param indexParams build parameters
   * @return total number of nucleotides read.
   * @throws IOException if error
   */
  static long indexThenSearchShortReads(final NgsParams params, final NgsHashLoop shl, final MapStatistics statistics, final CreateParams indexParams) throws IOException {
    Diagnostic.developerLog("index params: " + indexParams.toString());
    final HashFunctionFactory hashFunctionFactory = params.maskParams().maskFactory((int) params.getMaxReadLength());
    final IndexSet indexes = new IndexSet(params, indexParams, hashFunctionFactory.numberWindows());
    if (indexes.size() > INDEX_USAGE_REPORTING_THRESHOLD) {
      Diagnostic.warning("Selected parameters produce " + indexes.size() + " indexes (this is high and could be slow to run).");
    }
    final ReadCallImplementation rci = new ReadCallImplementation(indexes);
    final TemplateCallImplementation tci = new TemplateCallImplementation(params, indexParams.size(), indexes, null);

    final NgsHashFunction hf = hashFunctionFactory.create(rci, tci);

    final long numberReads = params.buildFirstParams().numberSequences() + (params.paired() ? params.buildSecondParams().numberSequences() : 0);
    hf.setReadSequences(numberReads);

    final long totalLength = index(params, shl, indexParams, indexes, hf);
    final OutputFilter filter = params.outputParams().outFilter();
    try (OutputProcessor outProcessor = filter.makeProcessor(params, statistics)) {
      tci.setOutputProcessor(outProcessor);

      search(params, shl, indexes, tci, hf);
      outProcessor.finish();
    }
    return totalLength;
  }

  private static long index(NgsParams params, NgsHashLoop shl, CreateParams indexParams, IndexSet indexes, NgsHashFunction hf) throws IOException {
    Diagnostic.developerLog("index start");
    long totalLength = 0;
    for (int pass = 1; pass <= (indexParams.compressHashes() ? 2 : 1); pass++) {
      totalLength = 0; //only count for one pass
      if (params.paired()) {
        final boolean cgFlip = params.buildFirstParams().reader().getPrereadType() == PrereadType.CG && params.buildFirstParams().reader().minLength() == CgUtils.CG_RAW_READ_LENGTH;
        final long l1 = shl.readLoop(params.buildFirstParams(), hf, ReadEncoder.PAIRED_FIRST, false);
        final long l2 = shl.readLoop(params.buildSecondParams(), hf, ReadEncoder.PAIRED_SECOND, cgFlip);
        totalLength += l1 + l2;

      } else {
        totalLength += shl.readLoop(params.buildFirstParams(), hf, ReadEncoder.SINGLE_END, false);
      }
      indexes.freeze(params.numberThreads());
    }
    return totalLength;
  }

  /**
   * Runs a search on the supplied template and indexes
   * @param params search parameters
   * @param shl the hash loop implementation to use
   * @param indexes indexes build upon the reads you are trying to map
   * @param tci where should hits be reported
   * @param hf the hash function to use
   * @throws IOException when the output feels poorly
   */
  public static void search(NgsParams params, NgsHashLoop shl, IndexSet indexes, TemplateCall tci, NgsHashFunction hf) throws IOException {
    //open the query reader
    //    NgsHashLoopImpl.READ_DELAY.reset();
    final Integer numberThreads = params.numberThreads();
    final ISequenceParams searchParams = params.searchParams();
    final int multiplier;
    if (numberThreads == 1) {
      multiplier = 1;
    } else {
      multiplier = params.threadMultiplier();
    }
    shl.templateLoopMultiCore(searchParams, hf, numberThreads, multiplier);

    //NgsHashLoopImpl.READ_DELAY.log("query");
    tci.logStatistics();
    for (int i = 0; i < indexes.size(); i++) {
      Diagnostic.userLog("Index[" + i + "] search performance " + LS + indexes.get(i).perfString());
    }
  }
}
