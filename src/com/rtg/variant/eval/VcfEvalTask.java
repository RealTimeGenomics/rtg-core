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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.CharBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SdfId;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.sam.SamRangeUtils;
import com.rtg.util.Pair;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.VcfUtils;


/**
 * Works out if calls are consistent with a baseline or not, always produces an ROC curve
 */
public final class VcfEvalTask extends ParamsTask<VcfEvalParams, NoStatistics> {

  private static final String FN_FILE_NAME = "fn.vcf";
  private static final String FP_FILE_NAME = "fp.vcf";
  private static final String TP_FILE_NAME = "tp.vcf";
  private static final String TPBASE_FILE_NAME = "tp-baseline.vcf";

  /** Filename used for the full ROC curve */
  public static final String FULL_ROC_FILE = "weighted_roc.tsv";
  /** Filename used for the homozygous curve */
  public static final String HOMOZYGOUS_FILE = "homozygous_roc.tsv";
  /** Filename used for the heterozygous curve */
  public static final String HETEROZYGOUS_FILE = "heterozygous_roc.tsv";
  /** Filename used for the simple ROC curve */
  public static final String SIMPLE_FILE = "simple_roc.tsv";
  /** Filename used for the complex ROC curve */
  public static final String COMPLEX_FILE = "complex_roc.tsv";
  /** Filename used for the homozygous simple ROC curve */
  public static final String HOMOZYGOUS_SIMPLE_FILE = "homozygous_simple_roc.tsv";
  /** Filename used for the homozygous complex ROC curve */
  public static final String HOMOZYGOUS_COMPLEX_FILE = "homozygous_complex_roc.tsv";
  /** Filename used for the heterozygous simple ROC curve */
  public static final String HETEROZYGOUS_SIMPLE_FILE = "heterozygous_simple_roc.tsv";
  /** Filename used for the heterozygous complex ROC curve */
  public static final String HETEROZYGOUS_COMPLEX_FILE = "heterozygous_complex_roc.tsv";

  private static BufferedReader getReader(final File f) throws IOException {
    return new BufferedReader(new InputStreamReader(FileUtils.createInputStream(f, true)));
  }

  protected VcfEvalTask(VcfEvalParams params, OutputStream reportStream, NoStatistics stats) {
    super(params, reportStream, stats, null);
  }

  @Override
  protected void exec() throws IOException {
    evaluateCalls(mParams);
  }

  /**
   * @param params the evaluation parameters.
   * @throws IOException when an IO exception occurs
   */
  static void evaluateCalls(VcfEvalParams params) throws IOException {
    SdfUtils.validateHasNames(params.templateFile());
    final SequencesReader templateSequences = SequencesReaderFactory.createMemorySequencesReader(params.templateFile(), true, LongRange.NONE);
    SdfUtils.validateNoDuplicates(templateSequences, params.templateFile(), false);

    final File baseline = params.baselineFile();
    final File calls = params.callsFile();
    final File output = params.directory();

    final boolean zip = params.outputParams().isCompressed();
    try (final BufferedReader baselineReader = getReader(baseline);
        final BufferedReader callReader = getReader(calls)) {
      checkHeader(baselineReader, callReader, templateSequences.getSdfId());
    }

    final List<Pair<String, Integer>> nameOrdering = new ArrayList<>();
    for (long i = 0; i < templateSequences.names().length(); i++) {
      nameOrdering.add(new Pair<>(templateSequences.names().name(i), templateSequences.length(i)));
    }

    final VariantSet variants = getVariants(nameOrdering, params, templateSequences);
    if (!output.exists() && !output.mkdirs()) {
      throw new IOException("Unable to create directory \"" + output.getPath() + "\"");
    }

    final boolean outputTpBase = params.outputBaselineTp();

    final String zipExt = zip ? FileUtils.GZ_SUFFIX : "";
    final File tpFile = new File(output, TP_FILE_NAME + zipExt);
    final File fpFile = new File(output, FP_FILE_NAME + zipExt);
    final File fnFile = new File(output, FN_FILE_NAME + zipExt);
    final File tpBaseFile = outputTpBase ? new File(output, TPBASE_FILE_NAME + zipExt) : null;

    try (final OutputStream tp = FileUtils.createOutputStream(tpFile, zip, false);
        final OutputStream fp = FileUtils.createOutputStream(fpFile, zip, false);
        final OutputStream fn = FileUtils.createOutputStream(fnFile, zip, false);
        final OutputStream tpBase = outputTpBase ? FileUtils.createOutputStream(tpBaseFile, zip, false) : null) {
      evaluateCalls(params, templateSequences, output, variants, tp, fp, fn, tpBase);
    }
    if (zip) {
      VcfUtils.createVcfTabixIndex(tpFile);
      VcfUtils.createVcfTabixIndex(fpFile);
      VcfUtils.createVcfTabixIndex(fnFile);
      if (tpBaseFile != null) {
        VcfUtils.createVcfTabixIndex(tpBaseFile);
      }
    }
  }

  private static void evaluateCalls(VcfEvalParams params, SequencesReader templateSequences, File output, VariantSet variants, OutputStream tp, OutputStream fp, OutputStream fn, OutputStream tpBase) throws IOException {
    final boolean zip = params.outputParams().isCompressed();
    tp.write(variants.calledHeader().toString().getBytes());
    fp.write(variants.calledHeader().toString().getBytes());
    fn.write(variants.baseLineHeader().toString().getBytes());
    if (tpBase != null) {
      tpBase.write(variants.baseLineHeader().toString().getBytes());
    }

    final PrereadNamesInterface names = templateSequences.names();
    final Map<String, Long> nameMap = new HashMap<>();
    for (long i = 0; i < names.length(); i++) {
      nameMap.put(names.name(i), i);
    }
    final EvalSynchronizer sync = new EvalSynchronizer(variants, tp, fp, fn, tpBase, params.baselineFile(), params.callsFile(), params.sortOrder());
    sync.mRoc.addFilter(RocFilter.ALL, new File(output, FULL_ROC_FILE));
    sync.mRoc.addFilter(RocFilter.HETEROZYGOUS, new File(output, HETEROZYGOUS_FILE));
    sync.mRoc.addFilter(RocFilter.HOMOZYGOUS, new File(output, HOMOZYGOUS_FILE));
    if (params.rtgStats()) {
      sync.mRoc.addFilter(RocFilter.SIMPLE, new File(output, SIMPLE_FILE));
      sync.mRoc.addFilter(RocFilter.COMPLEX, new File(output, COMPLEX_FILE));
      sync.mRoc.addFilter(RocFilter.HETEROZYGOUS_SIMPLE, new File(output, HETEROZYGOUS_SIMPLE_FILE));
      sync.mRoc.addFilter(RocFilter.HETEROZYGOUS_COMPLEX, new File(output, HETEROZYGOUS_COMPLEX_FILE));
      sync.mRoc.addFilter(RocFilter.HOMOZYGOUS_SIMPLE, new File(output, HOMOZYGOUS_SIMPLE_FILE));
      sync.mRoc.addFilter(RocFilter.HOMOZYGOUS_COMPLEX, new File(output, HOMOZYGOUS_COMPLEX_FILE));
    }

    final SimpleThreadPool threadPool = new SimpleThreadPool(params.numberThreads(), "VcfEval", true);
    threadPool.enableBasicProgress(templateSequences.numberSequences());
    for (int i = 0; i < templateSequences.numberSequences(); i++) {
      threadPool.execute(new SequenceEvaluator(sync, nameMap, templateSequences.copy()));
    }

    threadPool.terminate();

    if (variants.getNumberOfSkippedBaselineVariants() > 0) {
      Diagnostic.warning("There were " + variants.getNumberOfSkippedBaselineVariants() + " baseline variants skipped due to being too long, overlapping or starting outside the expected reference sequence length.");
    }
    if (variants.getNumberOfSkippedCalledVariants() > 0) {
      Diagnostic.warning("There were " + variants.getNumberOfSkippedCalledVariants() + " called variants skipped due to being too long, overlapping or starting outside the expected reference sequence length.");
    }
    if (sync.mRoc.getNumberOfIgnoredVariants() > 0) {
      Diagnostic.warning("There were " + sync.mRoc.getNumberOfIgnoredVariants() + " variants not included in ROC data files due to missing or invalid scores.");
    }
    Diagnostic.developerLog("Writing ROC");
    sync.mRoc.writeRocs(sync.mTruePositives + sync.mFalseNegatives, zip);
    if (params.outputSlopeFiles()) {
      produceSlopeFiles(params.directory(), zip, params.rtgStats());
    }
    writePhasingInfo(sync, params.directory());

    sync.mRoc.writeSummary(new File(params.directory(), CommonFlags.SUMMARY_FILE), sync.mTruePositives, sync.mFalsePositives, sync.mFalseNegatives);
  }

  private static void writePhasingInfo(EvalSynchronizer sync, File outDir) throws IOException {
    final File phasingFile = new File(outDir, "phasing.txt");
    FileUtils.stringToFile("Correct phasings: " + sync.getCorrectPhasings() + StringUtils.LS + "Incorrect phasings: " + sync.getMisPhasings() + StringUtils.LS + "Unresolvable phasings: " + sync.getUnphasable() + StringUtils.LS, phasingFile);
  }

  private static void produceSlopeFiles(File outDir, boolean zip, boolean rtgStats) throws IOException {
    final String suffix = zip ? FileUtils.GZ_SUFFIX : "";
    final File fullFile = new File(outDir, VcfEvalTask.FULL_ROC_FILE + suffix);
    produceSlopeFiles(fullFile, new File(outDir, "weighted_slope.tsv" + suffix), zip);
    final File heteroFile = new File(outDir, VcfEvalTask.HETEROZYGOUS_FILE + suffix);
    produceSlopeFiles(heteroFile, new File(outDir, "heterozygous_slope.tsv" + suffix), zip);
    final File homoFile = new File(outDir, VcfEvalTask.HOMOZYGOUS_FILE + suffix);
    produceSlopeFiles(homoFile, new File(outDir, "homozygous_slope.tsv" + suffix), zip);
    if (rtgStats) {
      final File simpleFile = new File(outDir, VcfEvalTask.SIMPLE_FILE + suffix);
      produceSlopeFiles(simpleFile, new File(outDir, "simple_slope.tsv" + suffix), zip);
      final File complexFile = new File(outDir, VcfEvalTask.COMPLEX_FILE + suffix);
      produceSlopeFiles(complexFile, new File(outDir, "complex_slope.tsv" + suffix), zip);
      final File heteroSimpleFile = new File(outDir, VcfEvalTask.HETEROZYGOUS_SIMPLE_FILE + suffix);
      produceSlopeFiles(heteroSimpleFile, new File(outDir, "heterozygous_simple_slope.tsv" + suffix), zip);
      final File heteroComplexFile = new File(outDir, VcfEvalTask.HETEROZYGOUS_COMPLEX_FILE + suffix);
      produceSlopeFiles(heteroComplexFile, new File(outDir, "heterozygous_complex_slope.tsv" + suffix), zip);
      final File homoSimpleFile = new File(outDir, VcfEvalTask.HOMOZYGOUS_SIMPLE_FILE + suffix);
      produceSlopeFiles(homoSimpleFile, new File(outDir, "homozygous_simple_slope.tsv" + suffix), zip);
      final File homoComplexFile = new File(outDir, VcfEvalTask.HOMOZYGOUS_COMPLEX_FILE + suffix);
      produceSlopeFiles(homoComplexFile, new File(outDir, "homozygous_complex_slope.tsv" + suffix), zip);
    }
  }

  private static void produceSlopeFiles(File input, File output, boolean zip) throws IOException {
    if (input.exists() && input.length() > 0) {
      try (final PrintStream printOut = new PrintStream(FileUtils.createOutputStream(output, zip));
          final InputStream in = zip ? FileUtils.createGzipInputStream(input, false) : FileUtils.createFileInputStream(input, false)) {
        RocSlope.writeSlope(in, printOut);
      }
    }
  }


  /**
   * Builds a variant set of the best type for the supplied files
   * @param nameOrdering the ordering of the names for output / processing
   * @param params the parameters
   * @param templateSequencesReader template sequences
   * @return a VariantSet for the provided files
   * @throws IOException if IO is broken
   */
  static VariantSet getVariants(Collection<Pair<String, Integer>> nameOrdering, VcfEvalParams params, SequencesReader templateSequencesReader) throws IOException {
    final File calls = params.callsFile();
    final File baseline = params.baselineFile();
    final String sampleName = params.sampleName();
    final RocSortValueExtractor extractor = getRocSortValueExtractor(params);
    final ReferenceRanges<String> ranges;
    if (params.bedRegionsFile() != null) {
      Diagnostic.developerLog("Loading BED regions");
      ranges = SamRangeUtils.createBedReferenceRanges(params.bedRegionsFile());
    } else if (params.restriction() != null) {
      ranges = SamRangeUtils.createExplicitReferenceRange(params.restriction());
    } else {
      ranges = SamRangeUtils.createFullReferenceRanges(templateSequencesReader);
    }
    return new TabixVcfRecordSet(baseline, calls, ranges, nameOrdering, sampleName, extractor, !params.useAllRecords(), params.squashPloidy(), params.maxLength());
  }

  private static RocSortValueExtractor getRocSortValueExtractor(VcfEvalParams params) {
    final RocScoreField fieldType;
    final String fieldName;
    final String scoreField = params.scoreField();
    if (scoreField != null) {
      if (scoreField.contains("=")) {
        final int pIndex = scoreField.indexOf('=');
        final String fieldTypeName = scoreField.substring(0, pIndex).toUpperCase(Locale.getDefault());
        try {
          fieldType = RocScoreField.valueOf(fieldTypeName);
        } catch (IllegalArgumentException e) {
          throw new NoTalkbackSlimException("Unrecognized field type \"" + fieldTypeName + "\", must be one of " + Arrays.toString(RocScoreField.values()));
        }
        fieldName = scoreField.substring(pIndex + 1);
      } else if (scoreField.equals(VcfUtils.QUAL)) {
        fieldType = RocScoreField.QUAL;
        fieldName = "UNUSED";
      } else {
        fieldType = RocScoreField.FORMAT;
        fieldName = scoreField;
      }
    } else {
      fieldType = RocScoreField.FORMAT;
      fieldName = VcfUtils.FORMAT_GENOTYPE_QUALITY;
    }
    return fieldType.getExtractor(fieldName, params.sortOrder());
  }

  static void checkHeader(BufferedReader baselineReader, BufferedReader callsReader, SdfId templateSdfId) throws IOException {
    final ArrayList<String> baselineHeader = readHeader(baselineReader);
    final ArrayList<String> callsHeader = readHeader(callsReader);

    if (baselineHeader.size() == 0) {
      throw new NoTalkbackSlimException("No header found in baseline file");
    } else if (callsHeader.size() == 0) {
      throw new NoTalkbackSlimException("No header found in calls file");
    }

    final SdfId baselineTemplateSdfId = getSdfId(baselineHeader);
    final SdfId callsTemplateSdfId = getSdfId(callsHeader);

    if (!baselineTemplateSdfId.check(templateSdfId)) {
      Diagnostic.warning("Template ID mismatch, baseline variants were not created from the given template");
    }

    if (!callsTemplateSdfId.check(templateSdfId)) {
      Diagnostic.warning("Template ID mismatch, called variants were not created from the given template");
    }

    if (!baselineTemplateSdfId.check(callsTemplateSdfId)) {
      Diagnostic.warning("Template ID mismatch, baseline and called variants were created with different templates");
    }
  }

  static ArrayList<String> readHeader(BufferedReader reader) throws IOException {
    if (reader == null) {
      return null;
    }
    final ArrayList<String> header = new ArrayList<>();
    String line;
    while ("#".equals(peek(reader)) && (line = reader.readLine()) != null) {
      header.add(line);
    }
    return header;
  }

  static String peek(BufferedReader reader) throws IOException {
    reader.mark(1);
    final CharBuffer buff = CharBuffer.allocate(1);
    if (reader.read(buff) == 1) {
      buff.rewind();
      reader.reset();
      return buff.toString();
    }
    return "";
  }

  static SdfId getSdfId(ArrayList<String> header) {
    if (header != null) {
      for (final String s : header) {
        if (s.startsWith("##TEMPLATE-SDF-ID=")) { //NOTE: this is brittle and VCF specific
          final String[] split = s.split("=");
          if (split.length != 2) {
            throw new NoTalkbackSlimException("Invalid header line : " + s);
          }
          final SdfId sdfId;
          try {
           sdfId = new SdfId(split[1]);
          } catch (final NumberFormatException ex) {
            throw new NoTalkbackSlimException("Invalid header line : " + s);
          }
          return sdfId;
        }
      }
    }
    return new SdfId(0);
  }

}

