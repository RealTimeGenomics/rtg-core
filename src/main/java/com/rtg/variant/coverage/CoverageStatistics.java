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
package com.rtg.variant.coverage;

import static com.rtg.util.StringUtils.LS;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.apache.velocity.tools.generic.NumberTool;

import com.reeltwo.plot.Datum2D;
import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.renderer.GraphicsRenderer;
import com.reeltwo.plot.ui.ImageWriter;
import com.rtg.launcher.AbstractStatistics;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.report.VelocityReportUtils;
import com.rtg.util.DoubleMultiSet;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.LongMultiSet;
import com.rtg.util.MathUtils;
import com.rtg.util.MultiSet;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.RangeMeta;
import com.rtg.util.io.FileUtils;

/**
 * Statistics object for tracking and printing coverage statistics.
 */
public class CoverageStatistics extends AbstractStatistics {

  /** Container for a row of coverage level output. */
  public static final class LevelDatum {
    final int mCoverageLevel;
    final long mCount;
    final double mPercent;
    final double mCumulative;

    private LevelDatum(int coverageLevel, long count, double percent, double cumulative) {
      mCoverageLevel = coverageLevel;
      mCount = count;
      mPercent = percent;
      mCumulative = cumulative;
    }

    public int getCoverageLevel() {
      return mCoverageLevel;
    }

    public long getCount() {
      return mCount;
    }

    public double getPercent() {
      return mPercent;
    }

    public double getCumulative() {
      return mCumulative;
    }
  }

  private static final int BUCKET_SIZE = Integer.getInteger("rtg.coverage.bucketsize", 1);

  private static final int BREADTH_DP = 4; // This is a fraction, so several DP is appropriate
  private static final int COVERAGE_DP = GlobalFlags.getIntegerValue(CoreGlobalFlags.COVERAGE_DP);
  private static final double RESIDUE_CUMULATIVE_PERCENT = 0.01; // Cumulative coverage below this level need not be reported

  private static final String COVERAGE_PNG_NAME = "coverage.png";
  private static final String CUMULATIVE_COVERAGE_PNG_NAME = "cumulative_coverage.png";
  private static final String STATS_TSV_NAME = "stats.tsv";
  private static final String LEVELS_TSV_NAME = "levels.tsv";

  private final boolean mDisableHtmlReport;
  private final int mFoldPercent;

  // Contains distinct labels to output (e.g. sequence names / region names / gene names
  private final LinkedHashSet<String> mCoverageNames = new LinkedHashSet<>();
  private final DoubleMultiSet<String> mTotalCoveragePerName = new DoubleMultiSet<>();
  private final LongMultiSet<String> mTotalLengthPerName = new LongMultiSet<>();
  private final LongMultiSet<String> mCoveredLengthPerName = new LongMultiSet<>();

  private final Map<RangeMeta<String>, CoverageSequenceStatistics> mOriginalRangeStatisticsMap = new HashMap<>();
  private List<RangeMeta<String>> mCurrentRange = Collections.emptyList();

  /**
   * Number of template positions whose coverage is in a given bucket.
   * Bucket 0 corresponds to coverage of exactly zero, while bucket i corresponds
   * to coverage from <code>(i - 1) * BUCKET_SIZE</code> (exclusive) up to
   * <code>i * BUCKET_SIZE</code> (inclusive).
   */
  private long[] mHistogram;

  /** the total number of bases across all regions, with overlaps being flattened. */
  private long mTotalBases;

  /** the number of covered bases across all regions, with overlaps being flattened. */
  private long mTotalCovered;

  /** the total coverage across all regions, with overlaps being flattened. */
  private double mTotalCoverage;

  private CoverageBedWriter mCoverageWriter;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   * @param disableHtmlReport true if HTML reporting should be disabled
   * @param targetFoldPercent the target percentage of bases when computing fold penalty statistic.
   */
  public CoverageStatistics(File outputDirectory, boolean disableHtmlReport, int targetFoldPercent) {
    super(outputDirectory);
    mFoldPercent = targetFoldPercent;
    mDisableHtmlReport = disableHtmlReport;
  }

  @Override
  protected String getStatistics() {
    return getStatistics(true);
  }

  String getStatistics(boolean summary) {
    final StringBuilder sb = new StringBuilder();
    final TextTable table = new TextTable(2, 2, TextTable.Align.RIGHT);

    if (summary) {
      sb.append("Coverage per region:").append(LS);
      appendRow(table, "depth", "breadth", "covered", "size", "name");
    } else {
      sb.append("#depth\tbreadth\tcovered\tsize\tname").append(LS);
    }
    final double[] cov = new double[mCoverageNames.size()];
    int k = 0;
    for (final String name : mCoverageNames) {
      final long size = mTotalLengthPerName.get(name);
      final double coverage = mTotalCoveragePerName.get(name);
      final long baseCount = mCoveredLengthPerName.get(name);
      cov[k++] = coverage / size;

      if (!summary || mCoverageNames.size() < 100) {  // show up to 99 labels in summary, otherwise only output the total line.
        final double breadth = size == 0 ? 0 : baseCount / (double) size;
        final double depth = size == 0 ? 0 : coverage / size;
        if (summary) {
          appendRow(table, Utils.realFormat(depth, COVERAGE_DP), Utils.realFormat(breadth, BREADTH_DP), baseCount, size, name);
        } else {
          sb.append(Utils.realFormat(depth, COVERAGE_DP)).append('\t').append(Utils.realFormat(breadth, BREADTH_DP)).append('\t')
                  .append(baseCount).append('\t').append(size).append('\t')
                  .append(name).append(LS);
        }
      }
    }
    final double totalDepth = mTotalBases == 0 ? 0.0 : mTotalCoverage / (double) mTotalBases;
    final double totalBreadth = mTotalBases == 0 ? 0.0 : mTotalCovered / (double) mTotalBases;

    if (summary) {
      appendRow(table, Utils.realFormat(totalDepth, COVERAGE_DP), Utils.realFormat(totalBreadth, BREADTH_DP), mTotalCovered, mTotalBases,
                     /*Utils.realFormat(nonNTotalDepth, DP), Utils.realFormat(nonNTotalBreadth, DP), mNonNTotalCovered, mNonNTotalBaseCount,*/ "all regions");
      table.toString(sb);
      final Double foldPenalty = foldPenalty(mFoldPercent);
      final double median = MathUtils.median(cov);
      if (foldPenalty != null || !Double.isNaN(median)) {
        sb.append(LS);
      }
      if (foldPenalty != null) {
        sb.append(String.format("Fold-%d Penalty: ", mFoldPercent)).append(Utils.realFormat(foldPenalty, 2)).append(LS);
      }
      if (!Double.isNaN(median)) {
        sb.append("Median depth: ").append(Utils.realFormat(median, COVERAGE_DP)).append(LS);
      }
    } else {
      sb.append(Utils.realFormat(totalDepth, COVERAGE_DP)).append('\t').append(Utils.realFormat(totalBreadth, BREADTH_DP)).append('\t')
              .append(mTotalCovered).append('\t').append(mTotalBases).append('\t')
              .append("all regions").append(LS);
    }
    return sb.toString();
  }

  void writeGraphs(List<LevelDatum> data, HtmlReportHelper hrh, Map<String, Object> velocityMap) throws IOException {

    final Graph2D graph = new Graph2D();
    graph.setGrid(true);

    final PointPlot2D plot = new PointPlot2D();
    plot.setLines(true);
    plot.setPoints(false);
    plot.setLineWidth(2);

    final Datum2D[] cumulativeCovData = new Point2D[data.size()];
    final Datum2D[] covData = new Point2D[data.size()];

    for (int i = 0; i < data.size(); ++i) {
      final LevelDatum datum = data.get(i);
      cumulativeCovData[i] = new Point2D(datum.mCoverageLevel, (float) datum.mCumulative); //%age cumulative
      covData[i] = new Point2D(datum.mCoverageLevel, (float) datum.mPercent); //%age
    }

    plot.setData(cumulativeCovData);
    plot.setColor(2);
    graph.addPlot(plot);

    graph.setLabel(Graph2D.X, "Coverage depth");
    graph.setLabel(Graph2D.Y, "Covered (cumulative %)");

    final ImageWriter iw = new ImageWriter(new GraphicsRenderer(new Color[] {Color.RED, new Color(0, 192, 0), new Color(0, 128, 255)}));
    iw.toImage(ImageWriter.ImageFormat.PNG, new File(hrh.getResourcesDir(), CUMULATIVE_COVERAGE_PNG_NAME), graph, 640, 480, null);
    velocityMap.put("cumulativeCoveragePng", CUMULATIVE_COVERAGE_PNG_NAME);

    final Graph2D graph2 = new Graph2D();
    graph2.setGrid(true);

    plot.setData(covData);
    graph2.addPlot(plot);

    graph2.setLabel(Graph2D.X, "Coverage depth");
    graph2.setLabel(Graph2D.Y, "Covered (%)");

    iw.toImage(ImageWriter.ImageFormat.PNG, new File(hrh.getResourcesDir(), COVERAGE_PNG_NAME), graph2, 640, 480, null);
    velocityMap.put("coveragePng", COVERAGE_PNG_NAME);
  }

  void writeLevels(List<LevelDatum> data, OutputStream os) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("#coverage_level\tcount\tpercent\tcumulative");
    sb.append(LS);

    for (final LevelDatum datum : data) {
      sb.append(datum.mCoverageLevel);
      sb.append('\t');
      sb.append(datum.mCount);
      sb.append('\t');
      sb.append(Utils.realFormat(datum.mPercent, 2)); //%age
      sb.append('\t');
      sb.append(Utils.realFormat(datum.mCumulative, 2)); //%age cumulative
      sb.append(LS);
    }
    os.write(sb.toString().getBytes());
  }

  private int getEffectiveLength() {
    // Since some runs can have extremely high coverage points (e.g. >100000x)
    // we sometimes want to bin multiple coverage points together so that the
    // output is smaller and the percentages reported are useful for plotting
    // etc.  However, it can be the case that the coverage distribution
    // contains a long tail of low-coverage counts and we don't want those to
    // affect the amount of binning that occurs.  Therefore we compute here
    // the number of coverage points that accounts for the majority of the data.
    long sum = 0;
    for (int k = 0; k < mHistogram.length; ++k) {
      final double cumulativeP = 100.0 * (mTotalBases - sum) / mTotalBases;
      sum += mHistogram[k];
      if (cumulativeP < RESIDUE_CUMULATIVE_PERCENT) {
        return k;
      }
    }
    return mHistogram.length;
  }

  List<LevelDatum> coverageLevels() {
    final List<LevelDatum> data = new ArrayList<>();
    if (mHistogram != null) {
      final int binSize = 1 + getEffectiveLength() / 1000;
      long sum = 0;
      for (int binStart = 0; binStart < mHistogram.length; binStart += binSize) {
        final double cumulativeP = 100.0 * (mTotalBases - sum) / mTotalBases;
        long binSum = 0;
        for (int j = 0; j < binSize && binStart + j < mHistogram.length; ++j) {
          binSum += mHistogram[binStart + j];
        }
        final double percent = 100.0 * binSum / mTotalBases;
        sum += binSum;
        data.add(new LevelDatum(binStart * BUCKET_SIZE, binSum, percent, cumulativeP));
        if (cumulativeP < RESIDUE_CUMULATIVE_PERCENT) {
          break;
        }
      }
    }
    return data;
  }

  private void appendRow(TextTable table, Object... parts) {
    final String[] values = new String[parts.length];
    for (int i = 0; i < parts.length; ++i) {
      values[i] = parts[i].toString();
    }
    table.addRow(values);
  }

  private static class CoverageSequenceStatistics {
    double mTotalCoverage = 0.0;
    int mBasesWithCoverage = 0;
    int mBases = 0;

    private void update(double nonSmoothCov, int minimumCoverageForBreadth) {
      mTotalCoverage += nonSmoothCov;
      ++mBases;
      if (nonSmoothCov >= minimumCoverageForBreadth) {
        ++mBasesWithCoverage;
      }
    }
  }

  void setPerRegionCoverageWriter(CoverageBedWriter writer) {
    mCoverageWriter = writer;
  }

  void setRange(String sequenceName, RangeList.RangeView<String> range) throws IOException {
    final List<RangeMeta<String>> newList = range == null ? Collections.emptyList() : range.getEnclosingRanges();

    //flush previous range stats, for any that are not present in the new list
    for (final RangeMeta<String> origRange : mCurrentRange) {
      if (!newList.contains(origRange)) {
        final CoverageSequenceStatistics stats = mOriginalRangeStatisticsMap.get(origRange);
        final String name = origRange.getMeta();
        if (mCoverageWriter != null && stats.mBases > 0) {
          mCoverageWriter.setRegionLabel(name);
          mCoverageWriter.finalCoverageRegion(sequenceName, origRange.getStart(), origRange.getEnd(), stats.mTotalCoverage / stats.mBases);
        }
        mTotalCoveragePerName.add(name, stats.mTotalCoverage);
        mTotalLengthPerName.add(name, stats.mBases);
        mCoveredLengthPerName.add(name, stats.mBasesWithCoverage);
        mOriginalRangeStatisticsMap.remove(origRange);
      }
    }

    // Create new range stats for any that aren't already present
    newList.stream().filter(origRange -> !mOriginalRangeStatisticsMap.containsKey(origRange)).forEach(origRange -> {
      mCoverageNames.add(origRange.getMeta());
      mOriginalRangeStatisticsMap.put(origRange, new CoverageSequenceStatistics());
    });

    mCurrentRange = newList;
  }

  /**
   * Update the coverage histogram and per-region-name statistics
   * This must be called once and only once for each base within any region.
   * @param currCoverage the current coverage
   * @param templateIsUnknown true if template is unknown (false if no template was supplied - treat as all known)
   * @param minimumCoverageForBreadth only increment covered bases stats if coverage above this value
   */
  public void updateCoverageHistogram(double currCoverage, boolean templateIsUnknown, int minimumCoverageForBreadth) {
    if (!templateIsUnknown) {
      final int bucket = ((int) currCoverage + BUCKET_SIZE - 1) / BUCKET_SIZE;
      if (mHistogram == null) {
        mHistogram = new long[bucket + 1];
      } else if (bucket >= mHistogram.length) {
        mHistogram = Arrays.copyOf(mHistogram, bucket + 1);
      }
      mHistogram[bucket]++;

      ++mTotalBases;

      if (currCoverage > 0) {
        if (currCoverage >= minimumCoverageForBreadth) {
          ++mTotalCovered;
        }
        mTotalCoverage += currCoverage;
      }
      for (final RangeMeta<String> origRange : mCurrentRange) {
        mOriginalRangeStatisticsMap.get(origRange).update(currCoverage, minimumCoverageForBreadth);
      }
    }
  }

  /**
   * Calculates the fold N penalty if possible. a measure of the non-uniformity of sequence
   * coverage: the amount of additional sequencing that would be necessary to ensure that N %
   * of target bases (in non-zero coverage targets) are covered to the current mean target coverage.
   * <br>
   * The fold penalty might not be able to be calculated if the binning results in no coverage value
   * near enough to the 20th percentile
   * @return the fold penalty, or null if it could not be calculated.
   * @param targetPct the percentage of target bases to be covered to mean coverage
   */
  public Double foldPenalty(int targetPct) {
    return foldPenalty(targetPct, mTotalBases, mTotalCoverage, mHistogram);
  }

  protected static Double foldPenalty(double targetPct, long totalBases, double totalCoverage, long[] histogram) {
    double sum = 0.0;
    if (histogram != null) {
      double lastCumulativePct = 100.0;
      for (int i = 0; i < histogram.length; ++i) {
        final double cumulativePct = 100.0 * (totalBases - sum) / totalBases;
        if (cumulativePct <= targetPct) {
          assert lastCumulativePct > targetPct;
          final double frac = (targetPct - cumulativePct) / (lastCumulativePct - cumulativePct);
          final double pct20Cov = (i - frac) * BUCKET_SIZE;
          return totalCoverage / totalBases / pct20Cov;
        }
        sum += histogram[i];
        lastCumulativePct = cumulativePct;
      }
    }
    return null;
  }

  @Override
  public void generateReport() throws IOException {
    final HtmlReportHelper hrh = getReportHelper();

    try (OutputStream statsFile = FileUtils.createOutputStream(new File(hrh.getBaseDir(), STATS_TSV_NAME))) {
      statsFile.write(getStatistics(false).getBytes());
    }

    final List<LevelDatum> levelsData = coverageLevels();

    try (OutputStream levelsOs = FileUtils.createOutputStream(new File(hrh.getBaseDir(), LEVELS_TSV_NAME))) {
      writeLevels(levelsData, levelsOs);
    }

    if (!mDisableHtmlReport) {
      final String m = ImageWriter.isImageWritingEnabled();
      if (m != null) {
        Diagnostic.warning("Skipping HTML report output, host OS is not correctly configured: " + m);
      } else {
        generateHtmlReport(levelsData, hrh);
      }
    }
  }

  void generateHtmlReport(List<LevelDatum> levelsData, HtmlReportHelper hrh) throws IOException {
    final Map<String, Object> velocityMap = new HashMap<>();

    velocityMap.put("resourceDir", hrh.getResourcesDirName());
    velocityMap.put("statsTsv", STATS_TSV_NAME);
    velocityMap.put("levelsTsv", LEVELS_TSV_NAME);
    velocityMap.put("levelsData", levelsData);

    writeGraphs(levelsData, hrh, velocityMap);

    velocityMap.put("numberTool", new NumberTool());
    velocityMap.put("commandLine", CommandLine.getCommandLine());
    String report = VelocityReportUtils.processTemplate("coverage.vm", velocityMap);
    report = VelocityReportUtils.wrapDefaultTemplate(report, "Coverage", hrh);

    FileUtils.stringToFile(report, hrh.getReportFile());
  }
}
