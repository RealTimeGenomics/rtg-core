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

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.STRING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collections;

import com.rtg.bed.BedUtils;
import com.rtg.bed.BedWriter;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Environment;
import com.rtg.util.MathUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.cnv.preprocess.AddGc;
import com.rtg.variant.cnv.preprocess.Column;
import com.rtg.variant.cnv.preprocess.GcNormalize;
import com.rtg.variant.cnv.preprocess.NumericColumn;
import com.rtg.variant.cnv.preprocess.RegionDataset;
import com.rtg.variant.cnv.preprocess.SimpleJoin;
import com.rtg.variant.cnv.preprocess.StringColumn;
import com.rtg.variant.cnv.preprocess.WeightedMedianNormalize;

/**
 * Provide construction of a panel sample for CNV calling.
 */
public class CnvPonBuildCli extends AbstractCli {

  private static final String VERSION_STRING = "#Version " + Environment.getVersion();
  private static final String CNV_PON_OUTPUT_VERSION = "v1.1";

  private static final String LABEL_COLUMN_NAME = "label-column-name";

  private static final String EXPANDED = "Xexpanded";

  // Name of the output column containing normalized coverage profile
  static final String NORMALIZED_COVERAGE_COLUMN = "normalized-coverage";

  @Override
  public String moduleName() {
    return "cnvponbuild";
  }

  @Override
  public String description() {
    return "build a typical CNV normalization sample from a panel of coverage output";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Construct a normalized coverage sample from a panel of coverage outputs, for use during segmentation.");
    CommonFlagCategories.setCategories(mFlags);
    CommonFlags.initForce(mFlags);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    CommonFlags.initReferenceTemplate(mFlags, true);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "BED output file").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(SegmentCli.GCBINS_FLAG, Integer.class, INT, "number of bins when applying GC correction", 10).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(SegmentCli.COV_COLUMN_NAME, String.class, STRING, "name of the coverage column in input data", SegmentCli.DEFAULT_COLUMN_NAME).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(LABEL_COLUMN_NAME, String.class, STRING, "if set, include region labels using the named column from the input data").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(EXPANDED, "if set, include additional columns in the output file").setCategory(SENSITIVITY_TUNING);
    final Flag<File> covFlag = mFlags.registerRequired(File.class, FILE, "coverage BED file").setCategory(INPUT_OUTPUT);
    covFlag.setMaxCount(Integer.MAX_VALUE);
    mFlags.setValidator(flags -> flags.checkInRange(SegmentCli.GCBINS_FLAG, 0, Integer.MAX_VALUE)
      && CommonFlags.validateOutputFile(flags, FileUtils.getOutputFileName((File) flags.getValue(OUTPUT_FLAG), !flags.isSet(NO_GZIP), BedUtils.BED_SUFFIX))
    );
  }

  private NumericColumn normalize(final File coverageFile, RegionDataset typicalSample, final boolean gcCorrect, final int gcbins) throws IOException {
    Diagnostic.info("Normalizing and G+C correcting " + coverageFile);
    final String coverageColumnName = (String) mFlags.getValue(SegmentCli.COV_COLUMN_NAME);
    final RegionDataset coverageData = RegionDataset.readFromBed(coverageFile, Collections.singletonList(new NumericColumn(coverageColumnName)));
    if (typicalSample.size() != coverageData.size()) {
      throw new NoTalkbackSlimException("Number of regions in " + coverageFile + " does not match a previous input file");
    }
    final int covCol = coverageData.columnId(coverageColumnName);
    if (covCol == -1) {
      throw new NoTalkbackSlimException("Could not find column named " + coverageColumnName + " in " + coverageFile);
    }
    if (gcCorrect) {
      new SimpleJoin(typicalSample, "").process(coverageData); // Join in pre-computed GC content columns
      new GcNormalize(covCol, gcbins).process(coverageData);
    }
    new WeightedMedianNormalize(coverageData.columns() - 1).process(coverageData);
    return coverageData.asNumeric(coverageData.columns() - 1);
  }

  private void writeBedHeader(final BedWriter bw) throws IOException {
    bw.writeln(VERSION_STRING + ", CNV panel BED output " + CNV_PON_OUTPUT_VERSION);
    if (CommandLine.getCommandLine() != null) {
      bw.writeComment("CL\t" + CommandLine.getCommandLine());
    }
    bw.writeComment("RUN-ID\t" + CommandLine.getRunId());
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final String labelColumn = mFlags.isSet(LABEL_COLUMN_NAME) ? (String) mFlags.getValue(LABEL_COLUMN_NAME) : null;
    try (final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG))) {
      final int gcbins = (Integer) mFlags.getValue(SegmentCli.GCBINS_FLAG);
      final AddGc gcCorrector = gcbins > 1 ? new AddGc(sr) : null;
      final RegionDataset typicalSample = RegionDataset.readFromBed((File) mFlags.getAnonymousValue(0), labelColumn == null ? Collections.emptyList() : Collections.singletonList(new StringColumn(labelColumn)));
      typicalSample.getColumns().removeIf((Column col) -> !col.getName().equals(labelColumn));
      final int cleanFirst = typicalSample.columns();
      if (gcCorrector != null) {
        Diagnostic.info("Computing per-region G+C content");
        gcCorrector.process(typicalSample);
      }

      final int first = typicalSample.columns();
      int sampleNum = 1;
      for (final Object coverageFile : mFlags.getAnonymousValues(0)) {
        final NumericColumn covData = normalize((File) coverageFile, typicalSample, gcCorrector != null, gcbins);
        covData.setName("wmednorm_" + sampleNum++);
        typicalSample.addColumn(covData);
      }
      final int last = typicalSample.columns();

      addPonMean(typicalSample, first, last);
      if (mFlags.isSet(EXPANDED)) {
        addPonMed(typicalSample, first, last);
      } else {
        for (int i = 0; i < (last - cleanFirst); i++) {
          typicalSample.getColumns().remove(cleanFirst);
        }
      }

      final boolean gzip = !mFlags.isSet(NO_GZIP);
      final File bedFile = FileUtils.getOutputFileName((File) mFlags.getValue(OUTPUT_FLAG), gzip, BedUtils.BED_SUFFIX);
      try (final BedWriter bw = new BedWriter(FileUtils.createOutputStream(bedFile))) {
        writeBedHeader(bw);
        typicalSample.write(bw);
      }
      final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
      if (gzip && index) {
        BedUtils.createBedTabixIndex(bedFile);
      }
    }
    return 0;
  }

  // Compute mean value for each ROI
  private void addPonMean(RegionDataset typicalSample, int first, int last) {
    final double[] sum = new double[typicalSample.size()];
    for (int i = first; i < last; i++) {
      final NumericColumn covData = typicalSample.asNumeric(i);
      for (int k = 0; k < sum.length; ++k) {
        sum[k] += covData.get(k);
      }
    }

    final int n = last - first;
    final NumericColumn col = new NumericColumn(NORMALIZED_COVERAGE_COLUMN);
    for (final double v : sum) {
      col.add(v / n);
    }
    typicalSample.addColumn(col);
  }

  // Compute median value for each ROI
  private void addPonMed(RegionDataset typicalSample, int first, int last) {
    final NumericColumn medcol = new NumericColumn("med-" + NORMALIZED_COVERAGE_COLUMN);
    final NumericColumn iqrcol = new NumericColumn("iqr-" + NORMALIZED_COVERAGE_COLUMN);
    final NumericColumn iqrnormcol = new NumericColumn("norm-iqr-" + NORMALIZED_COVERAGE_COLUMN);
    final double[] roi = new double[last - first];
    for (int k = 0; k < typicalSample.size(); ++k) {
      for (int i = first; i < last; i++) {
        roi[i - first] = typicalSample.asNumeric(i).get(k);
      }
      final double[] dist = MathUtils.quartiles(roi);
      final double iqr = dist[2] - dist[0];
      final double median = dist[1];
      medcol.add(median);
      iqrcol.add(iqr);
      iqrnormcol.add(median == 0 ? Double.POSITIVE_INFINITY : iqr / median);
    }
    typicalSample.addColumn(medcol);
    typicalSample.addColumn(iqrcol);
    typicalSample.addColumn(iqrnormcol);
  }
}
