/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.Properties;
import java.util.TreeSet;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedRangeLoader;
import com.rtg.bed.BedUtils;
import com.rtg.bed.BedWriter;
import com.rtg.bed.NamedBedRangeLoader;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.MathUtils;
import com.rtg.util.MultiSet;
import com.rtg.util.TextTable;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;
import com.rtg.variant.cnv.CnaType;
import com.rtg.variant.cnv.preprocess.AddGc;
import com.rtg.variant.cnv.preprocess.AddLog;
import com.rtg.variant.cnv.preprocess.AddRatio;
import com.rtg.variant.cnv.preprocess.Column;
import com.rtg.variant.cnv.preprocess.GcNormalize;
import com.rtg.variant.cnv.preprocess.IntersectJoin;
import com.rtg.variant.cnv.preprocess.NumericColumn;
import com.rtg.variant.cnv.preprocess.RegionColumn;
import com.rtg.variant.cnv.preprocess.RegionDataset;
import com.rtg.variant.cnv.preprocess.SimpleJoin;
import com.rtg.variant.cnv.preprocess.WeightedMedianNormalize;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.VcfWriterFactory;

/**
 * CNV segmentation task
 */
@TestClass("com.rtg.variant.cnv.segment.SegmentCliTest")
public class SegmentTask extends ParamsTask<SegmentParams, NoStatistics> {

  private final MultiSet<CnaType> mStatusCounts = new MultiSet<>();
  private SegmentVcfOutputFormatter mFormatter = null;
  private SequencesReader mReference = null;
  private RegionDataset mDataset;
  private int mDataCol;
  private CnvSummaryReport mReporter = null;
  private int mCaseCoverageCol = -1;
  private int mControlCoverageCol = -1;

  /**
   * @param params segmentation run parameters
   * @param reportStream stream to write statistics to.
   * @param stats statistics object to populate.
   */
  public SegmentTask(SegmentParams params, OutputStream reportStream, NoStatistics stats) {
    super(params, reportStream, stats, null);
  }

  @Override
  protected void exec() throws IOException {
    mReference = mParams.genome().reader();

    if (mParams.summaryRegionsFile() != null) {
      final ReferenceRanges<String> reportRegions = BedRangeLoader.getReferenceRanges(new NamedBedRangeLoader(), mParams.summaryRegionsFile());
      mReporter = new CnvSummaryReport(reportRegions);
    }

    if (mParams.precomputedColumn() != -1) {
      final File caseFile = mParams.caseFile();
      mDataset = RegionDataset.readFromBed(caseFile);
      mDataCol = mParams.precomputedColumn();
      Diagnostic.userLog("Using pre-computed column " + mDataCol + " for segmentation");

    } else if (mParams.panelFile() != null) {
      computeCasePanelDataset();
    } else {
      computeCaseControlDataset();
    }

    writeDataset();
    runSegmentation();
  }

  static final boolean LEGACY = false; // Simulate old behaviour where GcNormalization effectively happened after intersecting above-coverage rows

  // Input datasets are both coverage output, construct data based on (log) ratio with control
  private void computeCaseControlDataset() throws IOException {
    final double minCaseCoverage = mParams.minCaseCoverage();
    final int gcbins = mParams.gcBins();
    final String coverageColumnName = mParams.coverageColumnName();

    Diagnostic.userLog("Loading case");
    final File caseFile = mParams.caseFile();
    RegionDataset caseData = RegionDataset.readFromBed(caseFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (caseData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + caseFile);
    }
    mCaseCoverageCol = caseData.columnId(coverageColumnName);
    caseData.column(mCaseCoverageCol).setName("case_cover_raw");

    Diagnostic.userLog("Loading control");
    final File controlFile = mParams.controlFile();
    RegionDataset controlData = RegionDataset.readFromBed(controlFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (controlData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + controlData);
    }
    controlData.getColumns().removeIf((Column col) -> !col.getName().equals(coverageColumnName));
    mControlCoverageCol = controlData.columns() - 1;
    controlData.column(mControlCoverageCol).setName("ctrl_cover_raw");

    final RegionDataset gcDataset = new RegionDataset(controlData.regions());
    if (gcbins > 1) {
      Diagnostic.userLog("Computing per-region GC values");
      new AddGc(mReference).process(gcDataset);
    }

    if (mParams.minControlCoverage() != null) {
      // Apply minimum absolute control coverage filter
      final double minCtrlCoverage = mParams.minControlCoverage();
      final NumericColumn cc1 = controlData.asNumeric(mControlCoverageCol);
      controlData = controlData.filter(row -> cc1.get(row) >= minCtrlCoverage);
      Diagnostic.userLog("Filtered with minimum control coverage " + minCtrlCoverage + ", control dataset has " + controlData.size() + " rows");
    }

    final NumericColumn cc2 = caseData.asNumeric(mCaseCoverageCol);
    caseData = caseData.filter(row -> cc2.get(row) >= minCaseCoverage);
    Diagnostic.userLog("Filtered with minimum case coverage " + minCaseCoverage + ", case dataset has " + caseData.size() + " rows");

    if (gcbins > 1) {
      if (LEGACY) {
        new IntersectJoin(new RegionDataset(caseData.regions()), "", true, false).process(controlData);
        new IntersectJoin(new RegionDataset(controlData.regions()), "", true, false).process(caseData);
      }

      Diagnostic.userLog("Applying GC correction using " + gcbins + " bins");
      new IntersectJoin(gcDataset, "", true, false).process(caseData);
      mCaseCoverageCol = caseData.columns() - 1;
      new GcNormalize(mCaseCoverageCol, gcbins, "case_cover_gcnorm").process(caseData);
      mCaseCoverageCol = caseData.columns() - 1;


      new IntersectJoin(gcDataset, "", true, false).process(controlData);
      mControlCoverageCol = controlData.columns() - 1;
      new GcNormalize(mControlCoverageCol, gcbins, "ctrl_cover_gcnorm").process(controlData);
      mControlCoverageCol = controlData.columns() - 1;
    }

    Diagnostic.userLog("Applying weighted median normalization");
    new WeightedMedianNormalize(mCaseCoverageCol, "case_cover_wmednorm").process(caseData);
    mCaseCoverageCol = caseData.columns() - 1;
    new WeightedMedianNormalize(mControlCoverageCol, "ctrl_cover_wmednorm").process(controlData);
    mControlCoverageCol = controlData.columns() - 1;

    if (mParams.minControlCoverage() == null) {
      // Apply minimum normalized control coverage if appropriate
      final double minNormCtrlCoverage = mParams.minNormControlCoverage();
      final NumericColumn cc1 = controlData.asNumeric(mControlCoverageCol);
      controlData = controlData.filter(row -> cc1.get(row) >= minNormCtrlCoverage);
      Diagnostic.userLog("Filtered with minimum normalized control coverage " + minNormCtrlCoverage + ", control dataset has " + controlData.size() + " rows");
    }

    Diagnostic.userLog("Joining");
    new IntersectJoin(controlData, "").process(caseData);
    final RegionDataset filtered = caseData;
    mControlCoverageCol = filtered.columns() - 1;
    Diagnostic.userLog("Joined dataset has " + filtered.size() + " rows");

    Diagnostic.userLog("Computing ratio");
    //checkNonZero(filtered, caseCoverageCol);
    checkNonZero(filtered, mControlCoverageCol);
    new AddRatio(mCaseCoverageCol, mControlCoverageCol, "ratio_wmednorm").process(filtered);

    // Log
    Diagnostic.userLog("Computing log");
    new AddLog(filtered.columns() - 1, "ratio_wmednorm_log2").process(filtered);

    mDataset = filtered;
    mDataCol = mDataset.columns() - 1;
  }


  private void computeCasePanelDataset() throws IOException {
    final double minCaseCoverage = mParams.minCaseCoverage();
    final double minPanelCoverage = mParams.minNormControlCoverage();
    final int gcbins = mParams.gcBins();
    final String coverageColumnName = mParams.coverageColumnName();
    final String panelCoverageColumnName = mParams.panelCoverageColumnName();

    Diagnostic.userLog("Loading case");
    final File caseFile = mParams.caseFile();
    RegionDataset caseData = RegionDataset.readFromBed(caseFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (caseData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + caseFile);
    }

    mCaseCoverageCol = caseData.columnId(coverageColumnName);
    caseData.column(mCaseCoverageCol).setName("case_cover_raw");


    Diagnostic.userLog("Loading panel file");
    final File panelFile = mParams.panelFile();
    RegionDataset panelData = RegionDataset.readFromBed(panelFile, Collections.singletonList(new NumericColumn(panelCoverageColumnName)));
    if (panelData.columnId(panelCoverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + panelCoverageColumnName + " in " + panelFile);
    }
    panelData.getColumns().removeIf((Column col) -> !panelCoverageColumnName.equals(col.getName()));
    mControlCoverageCol = panelData.columns() - 1;
    panelData.column(mControlCoverageCol).setName("pon_cover_wmednorm");


    // Min coverage filters
    final NumericColumn cc1 = panelData.asNumeric(mControlCoverageCol);
    panelData = panelData.filter(row -> cc1.get(row) >= minPanelCoverage);
    Diagnostic.userLog("Filtered with minimum panel coverage " + minPanelCoverage + ", panel dataset has " + panelData.size() + " rows");

    final NumericColumn cc2 = caseData.asNumeric(mCaseCoverageCol);
    caseData = caseData.filter(row -> cc2.get(row) >= minCaseCoverage);
    Diagnostic.userLog("Filtered with minimum case coverage " + minCaseCoverage + ", case dataset has " + caseData.size() + " rows");

    if (gcbins > 1) {
      if (LEGACY) {
        new IntersectJoin(new RegionDataset(panelData.regions()), "", true, false).process(caseData);
      }

      final RegionDataset gcDataset = new RegionDataset(caseData.regions());
      Diagnostic.userLog("Computing per-region GC values");
      new AddGc(mReference).process(gcDataset);

      Diagnostic.userLog("Applying GC correction using " + gcbins + " bins");
      new SimpleJoin(gcDataset, "", true, false).process(caseData);
      mCaseCoverageCol = caseData.columns() - 1;
      new GcNormalize(mCaseCoverageCol, gcbins, "case_cover_gcnorm").process(caseData);
      mCaseCoverageCol = caseData.columns() - 1;
      // Panel is already gc normalized
    }

    Diagnostic.userLog("Applying weighted median normalization");
    new WeightedMedianNormalize(mCaseCoverageCol, "case_cover_wmednorm").process(caseData);
    mCaseCoverageCol = caseData.columns() - 1;
    // Panel is already median normalized


    Diagnostic.userLog("Joining");
    new IntersectJoin(panelData, "").process(caseData);
    final RegionDataset filtered = caseData;
    mControlCoverageCol = filtered.columns() - 1;
    Diagnostic.userLog("Joined dataset has " + filtered.size() + " rows");

    Diagnostic.userLog("Computing ratio");
    //checkNonZero(filtered, caseCoverageCol);
    checkNonZero(filtered, mControlCoverageCol);
    new AddRatio(mCaseCoverageCol, mControlCoverageCol, "ratio_wmednorm").process(filtered);

    // Log
    Diagnostic.userLog("Computing log");
    new AddLog(filtered.columns() - 1, "ratio_wmednorm_log2").process(filtered);

    mDataset = filtered;
    mDataCol = mDataset.columns() - 1;
  }

  private void writeDataset() throws IOException {
    final File bedFile = mParams.outputParams().outFile("unsegmented.bed");
    try (final BedWriter bw = new BedWriter(FileUtils.createOutputStream(bedFile))) {
      mDataset.write(bw);
    }
    if (mParams.outputParams().isCompressed()) {
      BedUtils.createBedTabixIndex(bedFile);
    }
  }


  private void runSegmentation() throws IOException {
    final int minBins = mParams.minBins();
    final SegmentScorer scorer = new EnergySegmentScorer(mParams.alpha(), mParams.aleph());
    final Collection<SegmentChain> sg = new ArrayList<>();
    final double refThreshold = mParams.minLogR();

    final File vcfFile = mParams.outFile("segments.vcf");

    final NumericColumn ratioCol = mDataset.asNumeric(mDataCol);
    final NumericColumn caseCol = mDataset.asNumeric(mCaseCoverageCol);
    final NumericColumn ctrlCol = mDataset.asNumeric(mControlCoverageCol);
    final RegionColumn regions = mDataset.regions();
    mFormatter = new SegmentVcfOutputFormatter(mReference, refThreshold, minBins, mParams.sampleName());
    try (final VcfWriter vw = new VcfWriterFactory().zip(mParams.outputParams().isCompressed()).make(mFormatter.header(), vcfFile)) {
      double prevMidPoint = -1;
      String prevSeqName = null;
      SegmentChain sc = new SegmentChain(scorer);

      for (int i = 0; i < mDataset.size(); ++i) {
        final SequenceNameLocus rec = regions.get(i);
        final String seqName = rec.getSequenceName();
        if (!seqName.equals(prevSeqName)) {
          prevMidPoint = -1;
          prevSeqName = seqName;
          sg.add(sc);
          sc = new SegmentChain(scorer);
        }
        final int start = rec.getStart();
        final int end = rec.getEnd();
        final double data = ratioCol.get(i);
        final long length = end - start;
        final double newMid = start + 0.5 * length;
        final double distPrev = prevMidPoint < 0 ? 0 : newMid - prevMidPoint;
        assert distPrev >= 0;
        prevMidPoint = newMid;
        sc.add(new Segment(seqName, start, end, data, distPrev, caseCol.get(i), ctrlCol.get(i)));
      }
      Diagnostic.progress("Starting segmentation");
      sg.add(sc);
      runSegmentation(vw, sg, mParams.beta());
    }

    if (mReporter != null) {
      Diagnostic.userLog("Writing region report");
      mReporter.report(vcfFile, mParams.outFile("summary.bed"));
    }

    writeSummary();
  }

  private void writeSummary() throws IOException {
    final TextTable summary = new TextTable(1, 0, TextTable.Align.RIGHT);

    summary.setAlignment(TextTable.Align.LEFT);
    summary.addRow("Total Segments:", String.valueOf(mStatusCounts.totalCount()));
    summary.addRow("Deletions:", String.valueOf(mStatusCounts.get(CnaType.DEL)));
    summary.addRow("Duplications:", String.valueOf(mStatusCounts.get(CnaType.DUP)));

    Diagnostic.userLog("SEGMENTATION SUMMARY");
    Diagnostic.userLog(summary.toString());
    try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(new File(mParams.directory(), CommonFlags.SUMMARY_FILE)), mReportStream))) {
      summaryOut.print(summary);
    }
  }

  private Collection<Segment> split(final Collection<SegmentChain> sg, final double deltaEnergyLimit) throws IOException {
    // Keep splitting until the energy limit or until the maximum number of segments,
    // whichever comes first.
    final int minSegments = mParams.minSegments();
    final int maxSegments = mParams.maxSegments();
    final Map<String, Long> seqNameMap = ReaderUtils.getSequenceNameMap(mReference);
    final Comparator<Segment> locusComparator = (a, b) -> {
      final int sc = seqNameMap.get(a.getSequenceName()).compareTo(seqNameMap.get(b.getSequenceName()));
      if (sc != 0) {
        return sc;
      }
      final int start = Integer.compare(a.getStart(), b.getStart());
      if (start != 0) {
        return start;
      }
      return Integer.compare(a.getEnd(), b.getEnd());
    };

    final TreeSet<Segment> orderByDeltaEnergyLimit = new TreeSet<>((a, b) -> {
      final int sc = Double.compare(b.deltaEnergy(), a.deltaEnergy());
      if (sc != 0) {
        return sc;
      }
      return locusComparator.compare(a, b);
    });
    for (final SegmentChain chain : sg) {
      orderByDeltaEnergyLimit.addAll(chain);
    }
    if (!orderByDeltaEnergyLimit.isEmpty()) {
      while (orderByDeltaEnergyLimit.size() < maxSegments) {
        final double dE = orderByDeltaEnergyLimit.first().deltaEnergy();
        if (dE < deltaEnergyLimit && orderByDeltaEnergyLimit.size() >= minSegments) {
          break; // Nothing further to be done, reached beta limit
        }
        // Split the top segment
        final Segment highest = orderByDeltaEnergyLimit.pollFirst();
        final Segment left = highest.left();
        final Segment right = highest.right();
        if (left == null || right == null) {
          break;
        }
        if (right.bins() == 1 && mParams.absorbSingletons()) {
          // Absorb right into left, pushing down the boundary change on the right
          final Segment newLeft = Segment.absorbRight(left, right);
          orderByDeltaEnergyLimit.add(newLeft); // right is discarded
        } else if (left.bins() == 1 && mParams.absorbSingletons()) {
          // Absorb left into right, pushing down the boundary change on the left
          final Segment newRight = Segment.absorbLeft(left, right);
          orderByDeltaEnergyLimit.add(newRight); // left is discarded
        } else {
          orderByDeltaEnergyLimit.add(left);
          orderByDeltaEnergyLimit.add(right);
        }
      }
    }

    // Reorder by genome locus in reference order
    final TreeSet<Segment> orderByReference = new TreeSet<>(locusComparator);
    orderByReference.addAll(orderByDeltaEnergyLimit);
    return orderByReference;
  }

  private void output(final VcfWriter writer, final Segment a, final Segment b, final Segment c) throws IOException {
    if (b != null) {
      final String currentName = b.getSequenceName();
      final Segment prev = a == null || !a.getSequenceName().equals(currentName) ? null : a;
      final Segment next = c == null || !c.getSequenceName().equals(currentName) ? null : c;
      final VcfRecord record = mFormatter.vcfRecord(prev, b, next);
      mStatusCounts.add(CnaType.valueOf(record));
      writer.write(record);
    }
  }

  private void runSegmentation(final VcfWriter writer, final Collection<SegmentChain> sg, final double beta) throws IOException {
    final double nu = nu(sg);
    final double sensitivityLimit = beta * nu;
    Diagnostic.userLog("Nu = " + Utils.realFormat(nu, 3));
    Diagnostic.userLog("Sensitivity limit = " + Utils.realFormat(sensitivityLimit, 3));

    for (final SegmentChain chain : sg) {
      chain.collapse();

      if (mParams.graphviz() && !chain.isEmpty()) {
        toGraphViz(chain, sensitivityLimit, mParams.file("chain-" + chain.get(0).getSequenceName() + ".dot"));
      }

    }
    final Collection<Segment> outputSegments = split(sg, sensitivityLimit);

    Segment b = null;
    Segment c = null;
    for (final Segment s : outputSegments) {
      final Segment a = b;
      b = c;
      c = s;
      output(writer, a, b, c);
    }
    output(writer, b, c, null);
  }

  private double nu(final Collection<SegmentChain> chains) {
    // total variability of data set
    double sum = 0;
    double sumSquares = 0;
    int size = 0;
    for (final SegmentChain chain : chains) {
      for (final Segment s : chain) {
        sum += s.sum();
        sumSquares += s.sumSquares();
      }
      size += chain.size();
    }
    final double mean = sum / size;
    final double var = sumSquares / size - mean * mean;
    return Math.sqrt(var);
  }

  private void checkNonZero(RegionDataset dataset, int col) {
    final NumericColumn c = dataset.asNumeric(col);
    for (int i = 0; i < c.size(); ++i) {
      if (MathUtils.approxEquals(c.get(i), 0.0, 0.000001)) {
        throw new NoTalkbackSlimException("Point with zero value in the data: " + dataset.regions().get(i));
      }
    }
  }


  private void toGraphViz(SegmentChain chain, double sl, File outfile) throws IOException {
    try (final LineWriter w = new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(outfile)))) {
      final StringBuilder sb = new StringBuilder();
      sb.append(GenomeRelationships.initGraph(new Properties(), "segmentation", "Segments for Chromosome " + chain.get(0).getSequenceName()));
      for (Segment s : chain) {
        graphVizSubtree(sb, s, sl);
      }
      sb.append("}\n");
      w.writeln(sb.toString());
    }
  }

  private void graphVizSubtree(StringBuilder sb, Segment parent, double sl) {
    sb.append(graphVizNodeId(parent)).append(' ').append(graphVizNode(parent, sl));
    graphVizSubtree(sb, parent, sl, parent.left());
    graphVizSubtree(sb, parent, sl, parent.right());
  }

  private void graphVizSubtree(StringBuilder sb, Segment parent, double sl, Segment child) {
    if (child != null) {
      sb.append(graphVizNodeId(parent)).append(" -> ").append(graphVizNodeId(child)).append(";\n");
      graphVizSubtree(sb, child, sl);
    }
  }

  private String graphVizNode(Segment s, double sl) {
    final String fillColor = s.bins() == 1 ? "white"
      : s.mean() > 1.1
      ? "dodgerblue"
      : s.mean() > 0.7
      ? "skyblue"
      : s.mean() < -0.8
      ? "pink"
      : s.mean() < -1.4
      ? "deeppink"
      : "grey";
    final String label = (s.getStart() + 1) + "\n"
      + s.getEnd() + "\n"
      + Utils.realFormat(s.mean(), 3) + " (" + s.bins() + ")";
    final String shape = "box";
    final String style = s.deltaEnergy() >= sl ? "filled" : "filled,rounded";
    return "["
      + "label=\""
      + label + "\""
      + ", shape=\"" + shape + "\""
      + ", style=\"" + style + "\", fillcolor=\"" + fillColor + "\""
      + "];\n";
  }

  private String graphVizNodeId(Segment s) {
    return (s.getStart() + 1) + "." + s.getEnd();
  }

}
