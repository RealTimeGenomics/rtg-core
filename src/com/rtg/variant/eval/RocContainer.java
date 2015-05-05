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
package com.rtg.variant.eval;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.rtg.ml.ContingencyTable;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 */
public class RocContainer {
  int mNoScoreVariants = 0;
  final List<SortedMap<Double, RocPoint>> mRocs = new ArrayList<>();
  final List<RocFilter> mFilters = new ArrayList<>();
  final List<File> mOutputFiles = new ArrayList<>();
  final Comparator<Double> mComparator;

  /**
   * Constructor
   * @param sortOrder the sort order that ensures that "good" scores are sorted before "bad" scores
   */
  public RocContainer(RocSortOrder sortOrder) {
    switch (sortOrder) {
      case ASCENDING:
        mComparator = new AscendingDoubleComparator();
        break;
      case DESCENDING:
      default:
        mComparator = new DescendingDoubleComparator();
        break;
    }
  }

  /**
   * @return the number of variants not included in ROC curves due to invalid score values
   */
  public int getNumberOfIgnoredVariants() {
    return mNoScoreVariants;
  }

  /**
   * Add an ROC output file
   * @param filter filter to apply to calls, or <code>null</code> for all calls
   * @param outputFile file to write output to.
   */
  public void addFilter(RocFilter filter, File outputFile) {
    mRocs.add(new TreeMap<Double, RocPoint>(mComparator));
    mFilters.add(filter);
    mOutputFiles.add(outputFile);
  }

  /**
   * add single result to ROC
   * @param primarySortValue normally the posterior score
   * @param weight weight of the call, 0.0 indicates a false positive with weight of 1
   * @param v detected variant
   */
  public void addRocLine(double primarySortValue, double weight, DetectedVariant v) {
    if (Double.isNaN(primarySortValue) || Double.isInfinite(primarySortValue)) {
      mNoScoreVariants++;
    } else {
      for (int i = 0; i < mFilters.size(); i++) {
        final SortedMap<Double, RocPoint> map = mRocs.get(i);
        final RocFilter filter = mFilters.get(i);
        if (v.filterAccept(filter)) {
          final RocPoint point = new RocPoint(primarySortValue, weight, weight > 0.0 ? 0 : 1);
          final RocPoint old = map.put(primarySortValue, point);
          if (old != null) {
            point.mTp += old.mTp;
            point.mFp += old.mFp;
          }
        }
      }
    }
  }

  private static final String HEADER = "#total baseline variants: ";
  private static final String HEADER2 = "#score true_positives false_positives".replaceAll(" ", "\t");
  private static void rocHeader(OutputStream out, int totalVariants) throws IOException {
    out.write((HEADER + totalVariants + StringUtils.LS).getBytes());
    out.write((HEADER2 + StringUtils.LS).getBytes());
  }

  /**
   * output ROC data to file
   * @param totalBaselineVariants total number of baseline variants
   * @param zip whether output should be compressed
   * @throws IOException if an IO error occurs
   */
  public void writeRocs(int totalBaselineVariants, boolean zip) throws IOException {
    for (int i = 0; i < mOutputFiles.size(); i++) {
      double tp = 0.0;
      int fp = 0;
      try (OutputStream os = FileUtils.createOutputStream(FileUtils.getZippedFileName(zip, mOutputFiles.get(i)), zip)) {
        rocHeader(os, totalBaselineVariants);
        for (final Map.Entry<Double, RocPoint> me : mRocs.get(i).entrySet()) {
          final RocPoint p = me.getValue();
          tp += p.mTp;
          fp += p.mFp;
          final StringBuilder sb = new StringBuilder();
          sb.append(Utils.realFormat(me.getKey(), 3));
          sb.append("\t");
          sb.append(Utils.realFormat(tp, 3));
          sb.append("\t");
          sb.append(Integer.toString(fp));
          sb.append(StringUtils.LS);
          os.write(sb.toString().getBytes());
        }
      }
    }
  }

  void writeSummary(File summaryFile, int truePositives, int falsePositives, int falseNegatives) throws IOException {
    final String summary;
    final int totalPositives = truePositives + falseNegatives;
    if (totalPositives > 0) {
      final TextTable table = new TextTable();
      table.addRow("Threshold", "True-pos", "False-pos", "False-neg", "Precision", "Sensitivity", "F-measure");
      table.addSeparator();

      final int totalRocPositives = truePositives + falseNegatives - mNoScoreVariants;
      final RocPoint best = (totalRocPositives > 0) ? bestFMeasure(totalRocPositives) : null;
      if (best == null) {
        Diagnostic.warning("Not enough ROC data to maximize F-measure, only un-thresholded statistics will be computed.");
      } else {
        addRow(table, Utils.realFormat(best.mThreshold, 3), best.mTp, best.mFp, totalRocPositives - best.mTp);
      }
      addRow(table, "None", truePositives, falsePositives, falseNegatives);
      summary = table.toString();
    } else {
      summary = "0 total baseline variants, no summary statistics available" + StringUtils.LS;
    }
    Diagnostic.info(summary);
    try (OutputStream os = FileUtils.createOutputStream(summaryFile, false)) {
      os.write(summary.getBytes());
    }
  }

  private static void addRow(final TextTable table, String threshold, double truePositive, double falsePositive, double falseNegative) {
    final double precision = ContingencyTable.precision(truePositive, falsePositive);
    final double recall = ContingencyTable.recall(truePositive, falseNegative);
    table.addRow(threshold,
      Integer.toString((int) truePositive),
      Integer.toString((int) falsePositive),
      Integer.toString((int) falseNegative),
      Utils.realFormat(precision, 4),
      Utils.realFormat(recall, 4),
      Utils.realFormat(ContingencyTable.fMeasure(precision, recall), 4));
  }

  // Find the threshold entry (from the ALL filter) that maximises f-measure, as an (arbitrary but) fair comparison point
  RocPoint bestFMeasure(int totalBaselineVariants) {
    for (int i = 0; i < mFilters.size(); i++) {
      final RocFilter filter = mFilters.get(i);
      if (filter == RocFilter.ALL) {
        double tp = 0.0;
        int fp = 0;
        double best = -1;
        RocPoint bestPoint = null;
        for (final Map.Entry<Double, RocPoint> me : mRocs.get(i).entrySet()) {
          final RocPoint p = me.getValue();
          tp += p.mTp;
          fp += p.mFp;
          final double precision = ContingencyTable.precision(tp, fp);
          final double recall = ContingencyTable.recall(tp, totalBaselineVariants - tp);
          final double fMeasure = ContingencyTable.fMeasure(precision, recall);
          if (fMeasure >= best) {
            best = fMeasure;
            bestPoint = new RocPoint(p.mThreshold, tp, fp);
          }
        }
        return bestPoint;
      }
    }
    return null;
  }

  private static final class RocPoint {
    double mThreshold;
    double mTp;
    int mFp;

    private RocPoint(double threshold, double tp, int fp) {
      mThreshold = threshold;
      mTp = tp;
      mFp = fp;
    }
  }

  private static class DescendingDoubleComparator implements Comparator<Double>, Serializable {
    @Override
    public int compare(Double o1, Double o2) {
      return o2.compareTo(o1);
    }
  }

  private static class AscendingDoubleComparator implements Comparator<Double>, Serializable {
    @Override
    public int compare(Double o1, Double o2) {
      return o1.compareTo(o2);
    }
  }
}
