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

package com.rtg.simulation;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfId;
import com.rtg.sam.SamFilterOptions;
import com.rtg.sam.SamUtils;
import com.rtg.util.Constants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.VcfReader;

import htsjdk.samtools.SamReader;

/**
 */
final class ReadSimEvalParams {

  private static final String FILE = "file";
  private static final String READ_SDF = "reads";
  private static final String VARIANCE = "variance";
  private static final String SCORE_HISTOGRAM = "score-histogram";
  private static final String MAPQ_HISTOGRAM = "mapq-histogram";
  private static final String MAPQ_ROC = "mapq-roc";

  private static final String MUTATIONS_VCF = "mutations-vcf";
  private static final String SAMPLE = "sample";

  private static final String VERBOSE = "verbose";

  // Debug type options
  private static final String MISMATCH = "Xmismatch";
  // private static final String GAP_EXTENTION = "gap-extention";
  private static final String GAP_OPENING = "Xgap-opening";
  private static final String DUMP_RECORDS = "Xdump-records";

  private final File mReadDir;
  private final File[] mSamFiles;
  private final File mMutationsVcf;

  private final String mSample;

  private final int mVariance;
  private final int mGapOpening;
  private final int mMismatch;


  private final boolean mIsPaired;
  private final boolean mVerbose;
  private final boolean mDumpRecords;
  private final boolean mScoreHistograms;
  private final boolean mMapQ;
  private final boolean mMapQRoc;

  static void initFlags(final CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('r', READ_SDF, File.class, CommonFlags.SDF, "SDF containing reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(flags);
    final Flag<File> inFlag = flags.registerRequired(File.class, FILE, "SAM/BAM format files");
    inFlag.setMinCount(1);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    inFlag.setCategory(CommonFlagCategories.INPUT_OUTPUT);
    CommonFlags.initNoMaxFile(flags);
    flags.registerOptional('M', MUTATIONS_VCF, File.class, FILE, "VCF file containing genomic mutations to be compensated for").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('S', SAMPLE, String.class, CommonFlags.STRING, "name of the sample to use from the mutation VCF file, will default to using the first sample in the file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('v', VARIANCE, Integer.class, CommonFlags.INT, "variation allowed in start position", 0).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('g', GAP_OPENING, Integer.class, CommonFlags.INT, "gap opening penalty", 2).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    //flags.registerOptional('e', GAP_EXTENTION, Integer.class, CommonFlags.INT, "gap extention penalty", (Integer) 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('m', MISMATCH, Integer.class, CommonFlags.INT, "mismatch penalty", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    SamFilterOptions.registerMinMapQFlag(flags);
    SamFilterOptions.registerMaxASMatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASUnmatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerExcludeMatedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);
    SamFilterOptions.registerExcludeDuplicatesFlag(flags);

    flags.registerOptional(DUMP_RECORDS, "dump annotated records to the screen").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(VERBOSE, "provide more detailed breakdown of stats").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(SCORE_HISTOGRAM, "output histogram of read alignment/generated scores").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MAPQ_HISTOGRAM, "output histogram of MAPQ scores").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MAPQ_ROC, "output ROC table with respect to MAPQ scores").setCategory(CommonFlagCategories.REPORTING);
    flags.setValidator(flags1 -> CommonFlags.validateOutputDirectory(flags1)
      && CommonFlags.checkFileList(flags1, null, null, Constants.MAX_OPEN_FILES)
      && CommonFlags.validateInputFile(flags1, MUTATIONS_VCF)
      && flags1.checkIf(SAMPLE, MUTATIONS_VCF)
      && flags1.checkInRange(VARIANCE, 0, Integer.MAX_VALUE));
  }


  ReadSimEvalParams(final CFlags flags) throws IOException {
    mReadDir = (File) flags.getValue(READ_SDF);
    mIsPaired = ReaderUtils.isPairedEndDirectory(mReadDir);
    mSamFiles = getSamFiles(flags);

    mMutationsVcf = (File) flags.getValue(MUTATIONS_VCF);
    String sample = (String) flags.getValue(SAMPLE);
    if (mMutationsVcf != null) {
      final VcfReader reader = VcfReader.openVcfReader(mMutationsVcf);
      reader.close();
      final List<String> samples = reader.getHeader().getSampleNames();
      if (sample == null) {
        sample = samples.get(0);
      } else if (!samples.contains(sample)) {
        throw new InvalidParamsException("Sample: \"" + sample + "\" not contained in mutations VCF file");
      }
    }
    mSample = sample;

    mVariance = (Integer) flags.getValue(VARIANCE);
    mGapOpening = (Integer) flags.getValue(GAP_OPENING);
    mMismatch = (Integer) flags.getValue(MISMATCH);
    mVerbose = flags.isSet(VERBOSE);
    mDumpRecords = flags.isSet(DUMP_RECORDS);
    mScoreHistograms = flags.isSet(SCORE_HISTOGRAM);
    mMapQ = flags.isSet(MAPQ_HISTOGRAM);
    mMapQRoc = flags.isSet(MAPQ_ROC);
    check(mReadDir, mSamFiles);
  }


  private static void check(File readDir, File[] samFiles) throws IOException {
    final SdfId sdfId = ReaderUtils.getSdfId(readDir);

    for (final File f : samFiles) {
      try (SamReader reader = SamUtils.makeSamReader(FileUtils.createInputStream(f, false))) {
        SamUtils.checkReadsGuid(reader.getFileHeader(), sdfId);
      }
    }
  }

  private File[] getSamFiles(CFlags flags) {
    final List<?> values = flags.getAnonymousValues(0);
    final List<File> files = new ArrayList<>(values.size());
    for (final Object f : values) {
      files.add((File) f);
    }
    return files.toArray(new File[0]);
  }

  File mutationsVcf() {
    return mMutationsVcf;
  }

  String sample() {
    return mSample;
  }

  int variance() {
    return mVariance;
  }

  int gapOpeningPenalty() {
    return mGapOpening;
  }

  int misMatchPenalty() {
    return mMismatch;
  }

  File readDirectory() {
    return mReadDir;
  }

  boolean isPaired() {
    return mIsPaired;
  }

  boolean verbose() {
    return mVerbose;
  }

  boolean dumpRecords() {
    return mDumpRecords;
  }

  boolean scoreHistograms() {
    return mScoreHistograms;
  }

  File[] samFiles() {
    return mSamFiles;
  }

  boolean mapQHistogram() {
    return mMapQ;
  }

  boolean mapQRoc() {
    return mMapQRoc;
  }

}
