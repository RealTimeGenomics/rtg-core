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
package com.rtg.variant.cnv.segment;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.FLOAT;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.SDF;
import static com.rtg.launcher.CommonFlags.STRING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.variant.cnv.segment.CnvPonBuildCli.NORMALIZED_COVERAGE_COLUMN;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collections;

import com.rtg.bed.BedUtils;
import com.rtg.bed.BedWriter;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.sam.SamRangeUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.MultiSet;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;
import com.rtg.variant.cnv.preprocess.AddGc;
import com.rtg.variant.cnv.preprocess.AddLog;
import com.rtg.variant.cnv.preprocess.AddRatio;
import com.rtg.variant.cnv.preprocess.Column;
import com.rtg.variant.cnv.preprocess.GcNormalize;
import com.rtg.variant.cnv.preprocess.NumericColumn;
import com.rtg.variant.cnv.preprocess.RegionColumn;
import com.rtg.variant.cnv.preprocess.RegionDataset;
import com.rtg.variant.cnv.preprocess.SimpleJoin;
import com.rtg.variant.cnv.preprocess.WeightedMedianNormalize;
import com.rtg.vcf.AsyncVcfWriter;
import com.rtg.vcf.DefaultVcfWriter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;

/**
 * Provide segmentation entry point.
 */
public class SegmentCli extends LoggedCli {

  private static final String ALEPH_FLAG = "aleph";
  private static final String ALPHA_FLAG = "alpha";
  private static final String BETA_FLAG = "beta";
  private static final String LIMIT_FLAG = "Xlimit";

  private static final String MIN_LOGR_FLAG = "min-logr";
  private static final String MIN_BINS_FLAG = "min-bins";

  private static final String CASE_FLAG = "case";
  private static final String CONTROL_FLAG = "control";
  private static final String PANEL_FLAG = "panel";
  private static final String SAMPLE_FLAG = "sample";

  private static final String SUMMARY_FLAG = "summary-regions";

  private static final String COLUMN_FLAG = "Xcolumn";
  static final String GCBINS_FLAG = "Xgcbins";
  private static final String MIN_CASE_COV_FLAG = "min-case-coverage";
  private static final String MIN_CTRL_COV_FLAG = "min-ctrl-coverage";
  static final String COV_COLUMN_NAME = "Xcolumn-name";

  static final String DEFAULT_COLUMN_NAME = "coverage";

  private final MultiSet<CnaType> mStatusCounts = new MultiSet<>();
  private SegmentVcfOutputFormatter mFormatter = null;
  private SequencesReader mReference = null;
  private VcfWriter mVcfOut = null;
  private RegionDataset mDataset;
  private int mDataCol;
  private CnvSummaryReport mReporter = null;


  @Override
  public String moduleName() {
    return "segment";
  }

  @Override
  public String description() {
    return "segment depth of coverage data to identify copy number alterations";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Segments depth of coverage data to identify regions of consistent copy number.");
    CommonFlagCategories.setCategories(mFlags);
    CommonFlags.initOutputDirFlag(mFlags);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    mFlags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, SDF, "SDF containing reference genome").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(CASE_FLAG, File.class, FILE, "BED file supplying per-region coverage data for the sample").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(CONTROL_FLAG, File.class, FILE, "BED file supplying per-region coverage data for control sample").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(PANEL_FLAG, File.class, FILE, "BED file supplying per-region panel normalized coverage data").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(SUMMARY_FLAG, File.class, FILE, "BED file supplying gene-scale regions to report CNV interactions with").setCategory(INPUT_OUTPUT);

    mFlags.registerOptional(ALEPH_FLAG, Double.class, FLOAT, "weighting factor for inter-segment distances during energy scoring", 0.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(ALPHA_FLAG, Double.class, FLOAT, "weighting factor for intra-segment distances during energy scoring", 0.001).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(BETA_FLAG, Double.class, FLOAT, "segmentation sensitivity factor", 0.5).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(LIMIT_FLAG, Integer.class, INT, "lower bound on the number of segments to be produced", 1).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_CASE_COV_FLAG, Double.class, FLOAT, "minimum case coverage required for a bin to be included in segmentation", 5.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_CTRL_COV_FLAG, Double.class, FLOAT, "minimum control coverage required for a bin to be included in segmentation", 300.0).setCategory(SENSITIVITY_TUNING);

    mFlags.registerOptional(MIN_BINS_FLAG, Integer.class, INT, "minimum number of bins required for copy number alteration to be called", 1).setCategory(REPORTING);
    mFlags.registerOptional('r', MIN_LOGR_FLAG, Double.class, FLOAT, "minimum (absolute) log ratio required for copy number alteration to be called", 0.2).setCategory(REPORTING);
    mFlags.registerOptional('s', SAMPLE_FLAG, String.class, STRING, "sample name to associate with the case sample", "SAMPLE").setCategory(REPORTING);

    mFlags.registerOptional('c', COLUMN_FLAG, Integer.class, INT, "just run segmentation directly on the specified data column from the case file, ignoring control data", 0).setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(GCBINS_FLAG, Integer.class, INT, "number of bins when applying GC correction", 10).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(COV_COLUMN_NAME, String.class, STRING, "name of the coverage column in input data", DEFAULT_COLUMN_NAME).setCategory(SENSITIVITY_TUNING);
    mFlags.setValidator(flags -> CommonFlags.validateOutputDirectory(flags)
      && flags.checkInRange(MIN_LOGR_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(COLUMN_FLAG, 0, Integer.MAX_VALUE)
      && flags.checkInRange(MIN_CASE_COV_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(MIN_CTRL_COV_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(GCBINS_FLAG, 0, Integer.MAX_VALUE)
      && CommonFlags.validateInputFile(flags, CASE_FLAG)
      && CommonFlags.validateInputFile(flags, CONTROL_FLAG)
      && CommonFlags.validateInputFile(flags, SUMMARY_FLAG)
      && flags.checkInRange(LIMIT_FLAG, 1, Integer.MAX_VALUE)
      && flags.checkXor(COLUMN_FLAG, CONTROL_FLAG, PANEL_FLAG)
    );
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }


  @Override
  protected int mainExec(OutputStream out, LogStream logStream) throws IOException {

    try (final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG))) {
      mReference = sr;

      if (mFlags.isSet(SUMMARY_FLAG)) {
        final ReferenceRanges<String> reportRegions = SamRangeUtils.createBedReferenceRanges((File) mFlags.getValue(SUMMARY_FLAG));
        mReporter = new CnvSummaryReport(reportRegions);
      }

      if (mFlags.isSet(COLUMN_FLAG)) {
        final File caseFile = (File) mFlags.getValue(CASE_FLAG);
        mDataset = RegionDataset.readFromBed(caseFile);
        mDataCol = (Integer) mFlags.getValue(COLUMN_FLAG);
        Diagnostic.userLog("Using pre-computed column " + mDataCol + " for segmentation");

      } else if (mFlags.isSet(PANEL_FLAG)) {
        computeCasePanelDataset();
      } else /*if (mFlags.isSet(CONTROL_FLAG))*/ {
        computeCaseControlDataset();
      }

      writeDataset();
      runSegmentation();
    }

    return 0;
  }

  // Input datasets are both coverage output, construct data based on (log) ratio with control
  private void computeCaseControlDataset() throws IOException {
    final double minCaseCoverage = (Double) mFlags.getValue(MIN_CASE_COV_FLAG);
    final double minCtrlCoverage = (Double) mFlags.getValue(MIN_CTRL_COV_FLAG);
    final int gcbins = (Integer) mFlags.getValue(GCBINS_FLAG);
    final String coverageColumnName = (String) mFlags.getValue(COV_COLUMN_NAME);

    Diagnostic.userLog("Loading case");
    final File caseFile = (File) mFlags.getValue(CASE_FLAG);
    final RegionDataset caseData = RegionDataset.readFromBed(caseFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (caseData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + caseFile);
    }
    //caseData.getColumns().removeIf((Column col) -> !col.getName().equals(COVERAGE_COLUMN_NAME));
    int caseCoverageCol = caseData.columnId(coverageColumnName);
    caseData.column(caseCoverageCol).setName("case_cover_raw");


    Diagnostic.userLog("Loading control");
    final File controlFile = (File) mFlags.getValue(CONTROL_FLAG);
    final RegionDataset controlData = RegionDataset.readFromBed(controlFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (controlData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + controlData);
    }
    controlData.getColumns().removeIf((Column col) -> !col.getName().equals(coverageColumnName));


    Diagnostic.userLog("Joining");
    new SimpleJoin(controlData, "control").process(caseData);
    RegionDataset filtered = caseData;
    int controlCoverageCol = filtered.columns() - 1;
    filtered.column(controlCoverageCol).setName("ctrl_cover_raw");


    // Min coverage filters
    final NumericColumn cc1 = filtered.asNumeric(controlCoverageCol);
    filtered = filtered.filter(row -> cc1.get(row) >= minCtrlCoverage);
    Diagnostic.userLog("Filtered with minimum control coverage " + minCtrlCoverage + ", dataset has " + filtered.size() + " rows");

    // TODO deal with deletions properly
    final NumericColumn cc2 = filtered.asNumeric(caseCoverageCol);
    filtered = filtered.filter(row -> cc2.get(row) >= minCaseCoverage);
    Diagnostic.userLog("Filtered with minimum case coverage " + minCaseCoverage + ", dataset has " + filtered.size() + " rows");

    if (gcbins > 1) {
      Diagnostic.userLog("Computing per-region GC values");
      new AddGc(mReference).process(filtered);
      Diagnostic.userLog("Applying GC correction using " + gcbins + " bins");
      new GcNormalize(caseCoverageCol, gcbins, "case_cover_gcnorm").process(filtered);
      caseCoverageCol = filtered.columns() - 1;
      new GcNormalize(controlCoverageCol, gcbins, "ctrl_cover_gcnorm").process(filtered);
      controlCoverageCol = filtered.columns() - 1;
    }

    Diagnostic.userLog("Applying weighted median normalization");
    new WeightedMedianNormalize(caseCoverageCol, "case_cover_wmednorm").process(filtered);
    caseCoverageCol = filtered.columns() - 1;
    new WeightedMedianNormalize(controlCoverageCol, "ctrl_cover_wmednorm").process(filtered);
    controlCoverageCol = filtered.columns() - 1;

    Diagnostic.userLog("Computing ratio");
    checkNonZero(filtered, caseCoverageCol);
    checkNonZero(filtered, controlCoverageCol);
    new AddRatio(caseCoverageCol, controlCoverageCol, "ratio_wmednorm").process(filtered);

    // Log
    Diagnostic.userLog("Computing log");
    new AddLog(filtered.columns() - 1, "ratio_wmednorm_log2").process(filtered);

    mDataset = filtered;
    mDataCol = mDataset.columns() - 1;
  }


  private void computeCasePanelDataset() throws IOException {
    final double minCaseCoverage = (Double) mFlags.getValue(MIN_CASE_COV_FLAG);
    final double minPanelCoverage = (Double) mFlags.getValue(MIN_CTRL_COV_FLAG); // Slight abuse of semantics
    final int gcbins = (Integer) mFlags.getValue(GCBINS_FLAG);
    final String coverageColumnName = (String) mFlags.getValue(COV_COLUMN_NAME);

    Diagnostic.userLog("Loading case");
    final File caseFile = (File) mFlags.getValue(CASE_FLAG);
    final RegionDataset caseData = RegionDataset.readFromBed(caseFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (caseData.columnId(coverageColumnName) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + caseFile);
    }
    //caseData.getColumns().removeIf((Column col) -> !col.getName().equals(COVERAGE_COLUMN_NAME));
    int caseCoverageCol = caseData.columnId(coverageColumnName);
    caseData.column(caseCoverageCol).setName("case_cover_raw");


    Diagnostic.userLog("Loading panel file");
    final File panelFile = (File) mFlags.getValue(PANEL_FLAG);
    final RegionDataset panelData = RegionDataset.readFromBed(panelFile, Collections.singletonList(new NumericColumn(NORMALIZED_COVERAGE_COLUMN)));
    if (panelData.columnId(NORMALIZED_COVERAGE_COLUMN) == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + NORMALIZED_COVERAGE_COLUMN + " in " + panelData);
    }
    panelData.getColumns().removeIf((Column col) -> !col.getName().equals(NORMALIZED_COVERAGE_COLUMN));


    Diagnostic.userLog("Joining");
    new SimpleJoin(panelData, "panel").process(caseData);
    RegionDataset filtered = caseData;
    final int panelCoverageCol = filtered.columns() - 1;
    filtered.column(panelCoverageCol).setName("pon_cover_wmednorm");
    Diagnostic.userLog("Joined dataset has " + filtered.size() + " rows");

    // Min coverage filters
    final NumericColumn cc1 = filtered.asNumeric(panelCoverageCol);
    filtered = filtered.filter(row -> cc1.get(row) >= minPanelCoverage);
    Diagnostic.userLog("Filtered with minimum panel coverage " + minPanelCoverage + ", dataset has " + filtered.size() + " rows");

    // TODO deal with deletions properly
    final NumericColumn cc2 = filtered.asNumeric(caseCoverageCol);
    filtered = filtered.filter(row -> cc2.get(row) >= minCaseCoverage);
    Diagnostic.userLog("Filtered with minimum case coverage " + minCaseCoverage + ", dataset has " + filtered.size() + " rows");

    if (gcbins > 1) {
      Diagnostic.userLog("Computing per-region GC values");
      new AddGc(mReference).process(filtered);
      Diagnostic.userLog("Applying GC correction using " + gcbins + " bins");
      new GcNormalize(caseCoverageCol, gcbins, "case_cover_gcnorm").process(filtered);
      caseCoverageCol = filtered.columns() - 1;
      // Panel is already gc normalized
    }

    Diagnostic.userLog("Applying weighted median normalization");
    new WeightedMedianNormalize(caseCoverageCol, "case_cover_wmednorm").process(filtered);
    caseCoverageCol = filtered.columns() - 1;
    // Panel is already median normalized

    Diagnostic.userLog("Computing ratio");
    checkNonZero(filtered, caseCoverageCol);
    checkNonZero(filtered, panelCoverageCol);
    new AddRatio(caseCoverageCol, panelCoverageCol, "ratio_wmednorm").process(filtered);

    // Log
    Diagnostic.userLog("Computing log");
    new AddLog(filtered.columns() - 1, "ratio_wmednorm_log2").process(filtered);

    mDataset = filtered;
    mDataCol = mDataset.columns() - 1;
  }

  private void writeDataset() throws IOException {
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
    final File bedFile = FileUtils.getZippedFileName(gzip, new File(outputDirectory(), "unsegmented.bed"));
    try (final BedWriter bw = new BedWriter(FileUtils.createOutputStream(bedFile, gzip))) {
      mDataset.write(bw);
    }
    if (gzip && index) {
      BedUtils.createBedTabixIndex(bedFile);
    }
  }


  private void runSegmentation() throws IOException {
    final int limit = (Integer) mFlags.getValue(LIMIT_FLAG);
    final int minBins = (Integer) mFlags.getValue(MIN_BINS_FLAG);
    final SegmentScorer scorer = new EnergySegmentScorer((Double) mFlags.getValue(ALPHA_FLAG), (Double) mFlags.getValue(ALEPH_FLAG));
    final SegmentChain sg = new SegmentChain(scorer, (Double) mFlags.getValue(BETA_FLAG));
    final double refThreshold = (Double) mFlags.getValue(MIN_LOGR_FLAG);

    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);

    final File vcfFile = VcfUtils.getZippedVcfFileName(gzip, new File(outputDirectory(), "segments.vcf"));

    final NumericColumn c = mDataset.asNumeric(mDataCol);
    final RegionColumn regions = mDataset.regions();
    mFormatter = new SegmentVcfOutputFormatter(mReference, refThreshold, minBins, (String) mFlags.getValue(SAMPLE_FLAG));
    try (final VcfWriter vw = new AsyncVcfWriter(new DefaultVcfWriter(mFormatter.header(), vcfFile, null, gzip, index))) {
      mVcfOut = vw;
      double prevMidPoint = -1;
      String prevSeq = null;

      for (int i = 0; i < mDataset.size(); ++i) {
        final SequenceNameLocus rec = regions.get(i);
        final String seq = rec.getSequenceName();
        final int start = rec.getStart();
        final int end = rec.getEnd();
        if (prevSeq != null && !seq.equals(prevSeq)) {
          Diagnostic.progress("Processing: " + prevSeq);
          // We are about to swap to a new sequence, so do all the processing for the records of prevSeq
          processSequence(prevSeq, sg, limit);
          prevMidPoint = -1;
          sg.clear();
        }
        prevSeq = seq;
        final double data = c.get(i);
        final long length = end - start;
        final double newMid = start + 0.5 * length;
        final double distPrev = prevMidPoint < 0 ? 0 : newMid - prevMidPoint;
        assert distPrev >= 0;
        prevMidPoint = newMid;
        sg.add(new Segment(start, end, data, distPrev));
      }
      if (prevSeq != null) {
        Diagnostic.progress("Processing: " + prevSeq);
        processSequence(prevSeq, sg, limit);
      }
    }

    if (mReporter != null) {
      Diagnostic.userLog("Writing region report");
      mReporter.report(vcfFile, FileUtils.getZippedFileName(gzip, new File(outputDirectory(), "summary.bed")));
    }

    writeSummary();
  }

  private void writeSummary() throws IOException {
    final TextTable summary = new TextTable(1, 0, TextTable.Align.RIGHT);

    summary.setAlignment(TextTable.Align.LEFT);
    summary.addRow("Total Segments:", "" + mStatusCounts.totalCount());
    summary.addRow("Deletions:", "" + mStatusCounts.get(CnaType.DEL));
    summary.addRow("Duplications:", "" + mStatusCounts.get(CnaType.DUP));

    Diagnostic.userLog("SEGMENTATION SUMMARY");
    Diagnostic.userLog(summary.toString());
    try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(new File(outputDirectory(), CommonFlags.SUMMARY_FILE), false), FileUtils.getStdoutAsOutputStream()))) {
      summaryOut.print(summary);
    }
  }

  private void processSequence(final String seqName, final SegmentChain sg, final int limit) throws IOException {
    sg.collapse(limit);
    for (int i = 0; i < sg.size(); ++i) {
      final Segment s = sg.get(i);
      final VcfRecord record = mFormatter.vcfRecord(seqName, s, i + 1 < sg.size() ? sg.get(i + 1) : null);
      mStatusCounts.add(CnaType.valueOf(record));
      mVcfOut.write(record);
    }
  }

  private void checkNonZero(RegionDataset dataset, int col) {
    final NumericColumn c = dataset.asNumeric(col);
    for (int i = 0; i < c.size(); ++i) {
      if (MathUtils.approxEquals(c.get(i), 0.0, 0.000001)) {
        throw new NoTalkbackSlimException("Point with zero value in the data: " + dataset.regions().get(i));
      }
    }
  }
}
