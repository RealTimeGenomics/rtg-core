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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.velocity.tools.generic.NumberTool;

import com.reeltwo.plot.Datum2D;
import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.renderer.GraphicsRenderer;
import com.reeltwo.plot.ui.ImageWriter;
import com.rtg.launcher.AbstractStatistics;
import com.rtg.report.VelocityReportUtils;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.MathUtils;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
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

  private final SortedSet<String> mCoverageRegions = new TreeSet<>();
  private final HashMap<String, Double> mTotalCoveragePerRegion = new HashMap<>();
  private final HashMap<String, Double> mTotalNonNCoveragePerRegion = new HashMap<>();
  private final HashMap<String, Integer> mRegionLengths = new HashMap<>();
  private final HashMap<String, Integer> mRegionCoveredBaseCounts = new HashMap<>();
  private final HashMap<String, Integer> mRegionNonNLengths = new HashMap<>();
  private final HashMap<String, Integer> mRegionNonNCoveredBaseCounts = new HashMap<>();

  private RangeList.RangeData<String> mCurrentRange = null;
  private Map<RangeList.RangeData<String>, CoverageSequenceStatistics> mOriginalRangeStatisticsMap = null;

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

  private Collection<String> mOutputRegions;

  //Decimal places for depth and breadth columns
  private static final int DP = 4;

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
    final Collection<String> sortOrder;
    if (mOutputRegions == null) {
      sortOrder = new ArrayList<>();
      sortOrder.addAll(mCoverageRegions);
    } else {
      sortOrder = mOutputRegions;
    }
    for (final String regionName : sortOrder) {
      final int size = mRegionLengths.get(regionName);
      final double coverage = mTotalCoveragePerRegion.get(regionName);
      final int baseCount = mRegionCoveredBaseCounts.get(regionName);
      final double breadth = baseCount / (double) size;
      final double depth = coverage / size;

      if (!summary || mCoverageRegions.size() < 100) {  // show up to 99 regions in summary, otherwise only output the total line.
        if (summary) {
          appendRow(table, Utils.realFormat(depth, DP), Utils.realFormat(breadth, DP), baseCount, size, regionName);
        } else {
          sb.append(Utils.realFormat(depth, DP)).append('\t').append(Utils.realFormat(breadth, DP)).append('\t')
                  .append(baseCount).append('\t').append(size).append('\t')
                  .append(regionName).append(LS);
        }
      }
    }
    final double totalDepth = mTotalBases == 0 ? 0.0 : mTotalCoverage / (double) mTotalBases;
    final double totalBreadth = mTotalBases == 0 ? 0.0 : mTotalCovered / (double) mTotalBases;

    if (summary) {
      appendRow(table, Utils.realFormat(totalDepth, DP), Utils.realFormat(totalBreadth, DP), mTotalCovered, mTotalBases,
                     /*Utils.realFormat(nonNTotalDepth, DP), Utils.realFormat(nonNTotalBreadth, DP), mNonNTotalCovered, mNonNTotalBaseCount,*/ "all regions");
      table.toString(sb);
    } else {
      sb.append(Utils.realFormat(totalDepth, DP)).append('\t').append(Utils.realFormat(totalBreadth, DP)).append('\t')
              .append(mTotalCovered).append('\t').append(mTotalBases).append('\t')
              .append("all regions").append(LS);
    }
    return sb.toString();
  }

  /**
   *  Set the region names to output, in the natural order of the supplied Collection
   *  @param  outputRegions the region names to output
   */
  void setOutputRegionNames(Collection<String> outputRegions) {
    mOutputRegions = outputRegions;
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

  /**
   * Add an average coverage value for a region to the statistics object.
   * @param regionName the name of the region
   * @param totalCoverage the total coverage
   * @param totalNonNCoverage the total coverage over non N bases
   * @param length the length of the region
   * @param nonNLength the non N length of the region
   */
  public void addAverageCoverage(String regionName, double totalCoverage, double totalNonNCoverage, int length, int nonNLength) {
    if (!mCoverageRegions.contains(regionName)) {
      mCoverageRegions.add(regionName);
    }
    incrementMap(mTotalCoveragePerRegion, regionName, totalCoverage);
    incrementMap(mTotalNonNCoveragePerRegion, regionName, totalNonNCoverage);
    incrementMap(mRegionLengths, regionName, length);
    incrementMap(mRegionNonNLengths, regionName, nonNLength);
  }

  void setRange(RangeList.RangeData<String> range) {
    if (mOriginalRangeStatisticsMap != null && mOriginalRangeStatisticsMap.size() > 0) {
      //flush previous range stats
      for (final RangeList.RangeData<String> origRange : mCurrentRange.getOriginalRanges()) {
        final String regionName = origRange.getMeta().get(0);
        final CoverageSequenceStatistics stats = mOriginalRangeStatisticsMap.get(origRange);
        addAverageCoverage(regionName, stats.mTotalNonNCoverage, stats.mTotalNonNCoverage, stats.mNonNLength, stats.mNonNLength);
        addCoveredBasesCount(regionName, stats.mNonNBasesWithCoverage);
        addNonNCoveredBasesCount(regionName, stats.mNonNBasesWithCoverage);
      }
      mOriginalRangeStatisticsMap.clear();
    }
    if (range != null) {
      if (mOriginalRangeStatisticsMap == null) {
        mOriginalRangeStatisticsMap = new HashMap<>();
      }
      mCurrentRange = range;
      for (final RangeList.RangeData<String> origRange : range.getOriginalRanges()) {
        mOriginalRangeStatisticsMap.put(origRange, new CoverageSequenceStatistics());
      }
    }
  }

  private static class CoverageSequenceStatistics {
    double mTotalNonNCoverage = 0.0;
    int mNonNBasesWithCoverage = 0;
    int mNonNLength = 0;

    private void update(double nonSmoothCov, int minimumCoverageForBreadth) {
      mTotalNonNCoverage += nonSmoothCov;
      mNonNLength++;
      if (nonSmoothCov >= minimumCoverageForBreadth) {
        mNonNBasesWithCoverage++;
      }
    }
  }

  private void incrementMap(Map<String, Double> map, String templateName, double newCount) {
    map.put(templateName, newCount + (map.containsKey(templateName) ? map.get(templateName) : 0));
  }
  private void incrementMap(Map<String, Integer> map, String templateName, int newCount) {
    map.put(templateName, newCount + (map.containsKey(templateName) ? map.get(templateName) : 0));
  }

  /**
   * Update the coverage histogram.
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
      for (final RangeList.RangeData<String> origRange : mCurrentRange.getOriginalRanges()) {
        mOriginalRangeStatisticsMap.get(origRange).update(currCoverage, minimumCoverageForBreadth);
      }
    }
  }

  /**
   * Add covered base count for a region
   * @param regionName name of the region
   * @param coveredBasesCount count of bases covered
   */
  public void addCoveredBasesCount(String regionName, Integer coveredBasesCount) {
    incrementMap(mRegionCoveredBaseCounts, regionName, coveredBasesCount);
  }

  /**
   * Add non N covered base count for a region
   * @param regionName name of the region
   * @param coveredBasesCount count of non N bases
   */
  public void addNonNCoveredBasesCount(String regionName, Integer coveredBasesCount) {
    incrementMap(mRegionNonNCoveredBaseCounts, regionName, coveredBasesCount);
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
