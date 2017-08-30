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
package com.rtg.variant.cnv;

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.variant.cnv.CnvProductParams.CnvProductParamsBuilder;

/**
 * Main class for CNV product
 */
public class CnvCli extends ParamsCli<CnvProductParams> {

  static final String BUCKET_SIZE_FLAG = "bucket-size";
  static final String INPUT_BASELINE_FLAG = "base-file";
  static final String INPUT_TEST_FLAG = "test-file";
  static final String INPUT_BASELINE_LIST_FLAG = "base-list-file";
  static final String INPUT_TEST_LIST_FLAG = "test-list-file";
  static final String ALLOW_MULTIPLE_ALIGNMENTS_PER_START_INDEX_FLAG = "Xallow-duplicate-start";
  static final String MAGIC_FLAG = "Xmagic-constant";
  static final String DIV_FACT_FLAG = "Xdivision-factor";
  static final String MUL_FACT_FLAG = "Xmultiplication-factor";
  static final String EXTRA_PENALTY_OFF = "Xextra-penalty-off";
  private static final String X_IGNORE_SAM_HEADER_INCOMPATIBILITY = "Xignore-incompatible-sam-headers";

  private static class CnvProductValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (flags.isSet(BUCKET_SIZE_FLAG) && (Integer) flags.getValue(BUCKET_SIZE_FLAG) < 1) {
        flags.setParseMessage("The bucket-size flag should be positive.");
        return false;
      }
      if (flags.isSet(CommonFlags.TEMPLATE_FLAG)) {
        final String genomePath = flags.getFlag(CommonFlags.TEMPLATE_FLAG).getValue().toString();
        final File genome = new File(genomePath);
        if (!genome.exists()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + genome.getPath() + "\", does not exist.");
          return false;
        }
        if (!genome.isDirectory()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + genome.getPath() + "\", is not an SDF.");
          return false;
        }
      }
      if (!flags.checkInRange(DIV_FACT_FLAG, 1.0, false, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!flags.checkInRange(MUL_FACT_FLAG, 1.0, false, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!flags.isSet(INPUT_BASELINE_FLAG) && !flags.isSet(INPUT_BASELINE_LIST_FLAG)) {
        flags.setParseMessage("Must set one of " + INPUT_BASELINE_LIST_FLAG + " or " + INPUT_BASELINE_FLAG);
        return false;
      }

      if (!flags.isSet(INPUT_TEST_FLAG) && !flags.isSet(INPUT_TEST_LIST_FLAG)) {
        flags.setParseMessage("Must set one of " + INPUT_TEST_LIST_FLAG + " or " + INPUT_TEST_FLAG);
        return false;
      }

      if (!CommonFlags.checkFileList(flags, INPUT_BASELINE_LIST_FLAG, INPUT_BASELINE_FLAG, Integer.MAX_VALUE)) {
        return false;
      }
      if (!CommonFlags.checkFileList(flags, INPUT_TEST_LIST_FLAG, INPUT_TEST_FLAG, Integer.MAX_VALUE)) {
        return false;
      }

      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }

      if (!SamFilterOptions.validateFilterFlags(flags, false)) {
        return false;
      }

      return true;
    }
  }

  @Override
  public String moduleName() {
    return "cnv";
  }

  @Override
  public String description() {
    return "call CNVs from paired SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Identifies copy number variation statistics and reports in a BED format file.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setValidator(new CnvProductValidator());
    mFlags.setName(applicationName() + " " + moduleName());
    final Flag<File> inFlag = mFlags.registerOptional('i', INPUT_BASELINE_FLAG, File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads for baseline");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> inFlag2 = mFlags.registerOptional('j', INPUT_TEST_FLAG, File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads for test");
    inFlag2.setCategory(INPUT_OUTPUT);
    inFlag2.setMinCount(0);
    inFlag2.setMaxCount(Integer.MAX_VALUE);
    CommonFlags.initNoMaxFile(mFlags);
    final Flag<File> listFlag1 = mFlags.registerOptional('I', INPUT_BASELINE_LIST_FLAG, File.class, CommonFlags.FILE, "file containing list of SAM/BAM format files (1 per line) containing mapped reads for baseline").setCategory(INPUT_OUTPUT);
    final Flag<File> listFlag2 = mFlags.registerOptional('J', INPUT_TEST_LIST_FLAG, File.class, CommonFlags.FILE, "file containing list of SAM/BAM format files (1 per line) containing mapped reads for test").setCategory(INPUT_OUTPUT);
    CommonFlags.initReferenceTemplate(mFlags, false);
    CommonFlags.initOutputDirFlag(mFlags);
    mFlags.registerOptional('b', BUCKET_SIZE_FLAG, Integer.class, CommonFlags.INT, "size of the buckets in the genome", 100).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(ALLOW_MULTIPLE_ALIGNMENTS_PER_START_INDEX_FLAG, "Count alignments starting at the same position individually").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MAGIC_FLAG, Double.class, CommonFlags.FLOAT, "Magic constant", 1.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(DIV_FACT_FLAG, Double.class, CommonFlags.FLOAT, "Division factor to use for calculating germline deletes", 3.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MUL_FACT_FLAG, Double.class, CommonFlags.FLOAT, "Multiplication factor to use for calculating germline deletes", 3.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(EXTRA_PENALTY_OFF, "Switch off extra penalty").setCategory(UTILITY);
    mFlags.registerOptional(X_IGNORE_SAM_HEADER_INCOMPATIBILITY, "ignore incompatible SAM headers when merging SAM results").setCategory(UTILITY);
    CommonFlags.initThreadsFlag(mFlags);
    CommonFlags.initNoGzip(mFlags);
    SamFilterOptions.registerMinMapQFlag(mFlags);
    SamFilterOptions.registerMaxHitsFlag(mFlags, 'c');
    SamFilterOptions.registerMaxASMatedFlag(mFlags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(mFlags, 'u');
    SamFilterOptions.registerExcludeMatedFlag(mFlags);
    SamFilterOptions.registerExcludeUnmatedFlag(mFlags);
    SamFilterOptions.registerRestrictionFlag(mFlags);
    mFlags.addRequiredSet(inFlag, inFlag2);
    mFlags.addRequiredSet(listFlag1, listFlag2);
  }

  @Override
  protected CnvProductParams makeParams() throws IOException {
    final CnvProductParamsBuilder builder = CnvProductParams.builder();
    builder.name(mFlags.getName());
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final OutputParams outParams = new OutputParams((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), mFlags.isSet(BuildCommon.PROGRESS_FLAG), gzip);
    builder.outputParams(outParams);
    final Collection<File> inputFilesBase = CommonFlags.getFileList(mFlags, INPUT_BASELINE_LIST_FLAG, INPUT_BASELINE_FLAG, false);
    Diagnostic.userLog("Base input SAM files: " + inputFilesBase);
    builder.mappedBase(inputFilesBase);
    final Collection<File> inputFilesTarget = CommonFlags.getFileList(mFlags, INPUT_TEST_LIST_FLAG, INPUT_TEST_FLAG, false);
    Diagnostic.userLog("Target input SAM files: " + inputFilesTarget);
    builder.mappedTarget(inputFilesTarget);
    builder.bucketSize((Integer) mFlags.getValue(BUCKET_SIZE_FLAG));
    builder.filterStartPositions(!mFlags.isSet(ALLOW_MULTIPLE_ALIGNMENTS_PER_START_INDEX_FLAG));
    builder.magicConstant((Double) mFlags.getValue(MAGIC_FLAG));
    builder.divisionFactor((Double) mFlags.getValue(DIV_FACT_FLAG));
    builder.multiplicationFactor((Double) mFlags.getValue(MUL_FACT_FLAG));
    builder.extraPenaltyOff(mFlags.isSet(EXTRA_PENALTY_OFF));
    builder.ignoreIncompatibleSamHeaders(mFlags.isSet(X_IGNORE_SAM_HEADER_INCOMPATIBILITY));
    if (mFlags.isSet(CommonFlags.TEMPLATE_FLAG)) {
      builder.genome(SequenceParams.builder().directory((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG)).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams());
    }
    builder.threads(CommonFlags.parseIOThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    final CnvProductParams localParams = builder.filterParams(SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnmapped(true).excludeUnplaced(true).create()).create();
    localParams.globalIntegrity();
    return localParams;

  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(CnvProductParams params, OutputStream out) {
    return new CnvProductTask(params, out);
  }

}
