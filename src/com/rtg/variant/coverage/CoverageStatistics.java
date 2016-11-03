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
import java.io.BufferedReader;
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
import com.rtg.report.VelocityReportUtils;
import com.rtg.util.DoubleMultiSet;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.MathUtils;
import com.rtg.util.MultiSet;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.io.FileUtils;

/**
 * Statistics object for tracking and printing coverage statistics.
 */
public class CoverageStatistics extends AbstractStatistics {

  private static final int BUCKET_SIZE = Integer.getInteger("rtg.coverage.bucketsize", 1);

  private static final String COVERAGE_PNG_NAME = "coverage.png";
  private static final String CUMULATIVE_COVERAGE_PNG_NAME = "cumulative_coverage.png";
  private static final String STATS_TSV_NAME = "stats.tsv";
  private static final String LEVELS_TSV_NAME = "levels.tsv";

  private final boolean mDisableHtmlReport;

  // Contains distinct labels to output (e.g. sequence names / region names / gene names
  private final LinkedHashSet<String> mCoverageNames = new LinkedHashSet<>();
  private final DoubleMultiSet<String> mTotalCoveragePerName = new DoubleMultiSet<>();
  private final MultiSet<String> mTotalLengthPerName = new MultiSet<>();
  private final MultiSet<String> mCoveredLengthPerName = new MultiSet<>();

  private final Map<RangeList.RangeData<String>, CoverageSequenceStatistics> mOriginalRangeStatisticsMap = new HashMap<>();
  private List<RangeList.RangeData<String>> mCurrentRange = Collections.emptyList();

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

  //Decimal places for depth and breadth columns
  private static final int DP = 4;
  private CoverageWriter mCoverageWriter;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   * @param disableHtmlReport true if HTML reporting should be disabled
   */
  public CoverageStatistics(File outputDirectory, boolean disableHtmlReport) {
    super(outputDirectory);
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
    for (final String name : mCoverageNames) {
      final int size = mTotalLengthPerName.get(name);
      final double coverage = mTotalCoveragePerName.get(name);
      final int baseCount = mCoveredLengthPerName.get(name);
      final double breadth = baseCount / (double) size;
      final double depth = coverage / size;

      if (!summary || mCoverageNames.size() < 100) {  // show up to 99 labels in summary, otherwise only output the total line.
        if (summary) {
          appendRow(table, Utils.realFormat(depth, DP), Utils.realFormat(breadth, DP), baseCount, size, name);
        } else {
          sb.append(Utils.realFormat(depth, DP)).append('\t').append(Utils.realFormat(breadth, DP)).append('\t')
                  .append(baseCount).append('\t').append(size).append('\t')
                  .append(name).append(LS);
        }
      }
    }
    final double totalDepth = mTotalBases == 0 ? 0.0 : mTotalCoverage / (double) mTotalBases;
    final double totalBreadth = mTotalBases == 0 ? 0.0 : mTotalCovered / (double) mTotalBases;

    if (summary) {
      appendRow(table, Utils.realFormat(totalDepth, DP), Utils.realFormat(totalBreadth, DP), mTotalCovered, mTotalBases,
                     /*Utils.realFormat(nonNTotalDepth, DP), Utils.realFormat(nonNTotalBreadth, DP), mNonNTotalCovered, mNonNTotalBaseCount,*/ "all regions");
      table.toString(sb);
      final Double fold80 = fold80();
      if (fold80 != null) {
        sb.append(LS);
        sb.append("Fold-80 Penalty: ").append(String.format("%.2f", fold80)).append(LS);
      }
    } else {
      sb.append(Utils.realFormat(totalDepth, DP)).append('\t').append(Utils.realFormat(totalBreadth, DP)).append('\t')
              .append(mTotalCovered).append('\t').append(mTotalBases).append('\t')
              .append("all regions").append(LS);
    }
    return sb.toString();
  }

  void writeGraphs(List<List<Double>> data, HtmlReportHelper hrh, Map<String, Object> velocityMap) throws IOException {

    final Graph2D graph = new Graph2D();
    graph.setGrid(true);

    final PointPlot2D plot = new PointPlot2D();
    plot.setLines(true);
    plot.setPoints(false);
    plot.setLineWidth(2);

    final Datum2D[] cumulativeCovData = new Point2D[data.size()];
    final Datum2D[] covData = new Point2D[data.size()];

    for (int i = 0; i < data.size(); i++) {
      final List<Double> datum = data.get(i);

      cumulativeCovData[i] = new Point2D(i, datum.get(2).floatValue()); //%age cumulative
      covData[i] = new Point2D(i, datum.get(1).floatValue()); //%age
    }

    plot.setData(cumulativeCovData);
    plot.setColor(2);
    graph.addPlot(plot);

    graph.setLabel(Graph2D.X, "Coverage depth");
    graph.setLabel(Graph2D.Y, "Covered (cumulative %)");

    final ImageWriter iw = new ImageWriter(new GraphicsRenderer(new Color[] {Color.RED, new Color(0, 192, 0), new Color(0, 128, 255)}));
    iw.toPNG(new File(hrh.getResourcesDir(), CUMULATIVE_COVERAGE_PNG_NAME), graph, 640, 480, null);
    velocityMap.put("cumulativeCoveragePng", CUMULATIVE_COVERAGE_PNG_NAME);

    final Graph2D graph2 = new Graph2D();
    graph2.setGrid(true);

    plot.setData(covData);
    graph2.addPlot(plot);

    graph2.setLabel(Graph2D.X, "Coverage depth");
    graph2.setLabel(Graph2D.Y, "Covered (%)");

    iw.toPNG(new File(hrh.getResourcesDir(), COVERAGE_PNG_NAME), graph2, 640, 480, null);
    velocityMap.put("coveragePng", COVERAGE_PNG_NAME);
  }

  void writeLevels(List<List<Double>> data, OutputStream os) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("#coverage_level\tcount\t%age\t%cumulative");
    sb.append(LS);

    String leader = "";
    for (int i = 0; i < data.size(); i++) {

      final List<Double> datum = data.get(i);
      sb.append(leader);
      sb.append(i * BUCKET_SIZE);
      sb.append('\t');
      sb.append(datum.get(0).longValue()); //count
      sb.append('\t');
      sb.append(Utils.realFormat(datum.get(1), 2)); //%age
      sb.append('\t');
      sb.append(Utils.realFormat(datum.get(2), 2)); //%age cumulative
      sb.append(LS);

      if (BUCKET_SIZE > 1) {
        leader = (i * BUCKET_SIZE + 1) + "..";
      }
    }
    os.write(sb.toString().getBytes());
  }

  List<List<Double>> coverageLevels() {
    final List<List<Double>> data = new ArrayList<>();

    if (mHistogram != null) {
      int lastNonZero = 0;
      final double[] nums = new double[mHistogram.length];
      for (int i = 0; i < mHistogram.length; i++) {
        final double percent = 100.0 * mHistogram[i] / mTotalBases;
//        final String num = Utils.realFormat(percent, 2);
        nums[i] = percent;
        if (!MathUtils.approxEquals(percent, 0.00, 0.01)) {
          lastNonZero = i;
        }
      }

      long sum = 0;
      for (int i = 0; i <= lastNonZero; i++) {
        final double cumulativeP = 100.0 * (mTotalBases - sum) / mTotalBases;
        sum += mHistogram[i];
        final List<Double> datum = new ArrayList<>();
        datum.add((double) mHistogram[i]);
        datum.add(nums[i]);
        datum.add(cumulativeP);
        data.add(datum);
      }
    }
    return data;
  }

  private void appendRow(TextTable table, Object... parts) {
    final String[] values = new String[parts.length];
    for (int i = 0; i < parts.length; i++) {
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
      mBases++;
      if (nonSmoothCov >= minimumCoverageForBreadth) {
        mBasesWithCoverage++;
      }
    }
  }

  void setPerRegionCoverageWriter(CoverageWriter writer) {
    mCoverageWriter = writer;
  }

  void setRange(String sequenceName, RangeList.RangeData<String> range) throws IOException {
    final List<RangeList.RangeData<String>> newList = range == null ? Collections.emptyList() : range.getOriginalRanges();

    //flush previous range stats, for any that are not present in the new list
    for (final RangeList.RangeData<String> origRange : mCurrentRange) {
      if (!newList.contains(origRange)) {
        final CoverageSequenceStatistics stats = mOriginalRangeStatisticsMap.get(origRange);
        final String name = origRange.getMeta().get(0);
        if (mCoverageWriter != null) {
          mCoverageWriter.setBedLabel(name);
          mCoverageWriter.finalCoverageRegion(sequenceName, origRange.getStart(), origRange.getEnd(), (int) MathUtils.round(stats.mTotalCoverage / stats.mBases));
        }
        mTotalCoveragePerName.add(name, stats.mTotalCoverage);
        mTotalLengthPerName.add(name, stats.mBases);
        mCoveredLengthPerName.add(name, stats.mBasesWithCoverage);
        mOriginalRangeStatisticsMap.remove(origRange);
      }
    }

    // Create new range stats for any that aren't already present
    newList.stream().filter(origRange -> !mOriginalRangeStatisticsMap.containsKey(origRange)).forEach(origRange -> {
      mCoverageNames.add(origRange.getMeta().get(0));
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

      mTotalBases++;

      if (currCoverage > 0) {
        if (currCoverage >= minimumCoverageForBreadth) {
          mTotalCovered++;
        }
        mTotalCoverage += currCoverage;
      }
      for (final RangeList.RangeData<String> origRange : mCurrentRange) {
        mOriginalRangeStatisticsMap.get(origRange).update(currCoverage, minimumCoverageForBreadth);
      }
    }
  }

  /**
   * Calculates the fold 80 penalty if possible. a measure of the non-uniformity of sequence
   * coverage: the amount of additional sequencing that would be necessary to ensure that 80%
   * of target bases (in non-zero coverage targets) are covered to the current mean target coverage.
   * <br>
   * The fold 80 penalty might not be able to be calculated if the binning results in no coverage value
   * near enough to the 20th percentile
   * @return the fold 80 penalty, or null if it could not be calculated.
   */
  public Double fold80() {
    double sum = 0.0;
    if (mHistogram != null) {
      for (int i = 0; i < mHistogram.length; i++) {
        final double cumulativeP = 100.0 * (mTotalBases - sum) / mTotalBases;
        sum += mHistogram[i];

        if (cumulativeP <= 80.0 + 1e-10) {
          if (MathUtils.approxEquals(cumulativeP, 80.0, 1e-1)) {
            final int pct20Cov = i * BUCKET_SIZE;
            return mTotalCoverage / mTotalBases / pct20Cov;
          }
          break;
        }
      }
    }
    return null;
  }

  @Override
  public void generateReport() throws IOException {
    final HtmlReportHelper hrh = getReportHelper();

    try (OutputStream statsFile = FileUtils.createOutputStream(new File(hrh.getBaseDir(), STATS_TSV_NAME), false)) {
      statsFile.write(getStatistics(false).getBytes());
    }

    final List<List<Double>> levelsData = coverageLevels();

    try (OutputStream levelsOs = FileUtils.createOutputStream(new File(hrh.getBaseDir(), LEVELS_TSV_NAME), false)) {
      writeLevels(levelsData, levelsOs);
    }

    if (!mDisableHtmlReport) {
      generateHtmlReport(levelsData, hrh);
    }
  }

  void generateHtmlReport(List<List<Double>> levelsData, HtmlReportHelper hrh) throws IOException {
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

  private List<List<Double>> readLevelsData(File levelsFile) throws IOException {
    final List<List<Double>> data = new ArrayList<>();
    try (final BufferedReader br = new BufferedReader(FileUtils.createReader(levelsFile, false))) {
      String line;
      while ((line = br.readLine()) != null) {
        if (!line.startsWith("#")) {
          final String[] split = StringUtils.split(line, '\t');
          if (split.length == 4) {
            final List<Double> datum = new ArrayList<>();
            datum.add(Double.parseDouble(split[1]));
            datum.add(Double.parseDouble(split[2]));
            datum.add(Double.parseDouble(split[3]));
            data.add(datum);
          }
        }
      }
    }
    return data;
  }

  /**
   * @param args arg arg arg
   * @throws Exception if bad
   */
  public static void main(String[] args) throws Exception {
    if (args.length != 1) {
      System.err.println("Usage: CoverageStatistics <coveragedir>");
      System.exit(1);
    }
    final File dir = new File(args[0]);
    if (!dir.exists() || !dir.isDirectory()) {
      throw new IOException("Could not find " + args[0] + " or was not a directory");
    }

    final File levels = new File(dir, LEVELS_TSV_NAME);
    if (!levels.exists()) {
      throw new IOException("Could not find file: " + args[0] + StringUtils.FS + LEVELS_TSV_NAME);
    }

    final CoverageStatistics stats = new CoverageStatistics(dir, false);

    //read stats levels file into double list list
    final List<List<Double>> data = stats.readLevelsData(levels);

    stats.generateHtmlReport(data, stats.getReportHelper());

  }
}
