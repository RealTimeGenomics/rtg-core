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
package com.rtg.sam;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import com.rtg.calibrate.SamCalibrationInputs;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;

import htsjdk.samtools.SAMFileHeader;

/**
 */
public class SamMergeCli extends AbstractCli {

  private static final String MODULE_NAME = "sammerge";

  private static final String REFERENCE_FLAG = "reference";
  private static final String LEGACY_CIGARS = "legacy-cigars";
  private static final String X_ALTERNATE_SAM_HEADER = "Xalternate-sam-header";

  private static class SamFileMergeValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final File o = (File) flags.getValue(CommonFlags.OUTPUT_FLAG);
      if (o != null) {
        final boolean isbam = o.getName().endsWith(SamUtils.BAM_SUFFIX);
        final File output = isbam ? o : SamUtils.getZippedSamFileName(!flags.isSet(CommonFlags.NO_GZIP), o);
        if (output.exists()) {
          flags.setParseMessage("The file \"" + output.getPath() + "\" already exists. Please remove it first or choose a different file");
          return false;
        }
      }
      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }

      if (!SamFilterOptions.validateFilterFlags(flags, true)) {
        return false;
      }
      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "merge sorted SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Merges and filters coordinate-sorted SAM/BAM files.");
    mFlags.registerOptional('t', CommonFlags.TEMPLATE_FLAG, File.class, "SDF", "SDF containing reference used to decode CRAM input").setCategory(INPUT_OUTPUT);
    CommonFlagCategories.setCategories(mFlags);
    final Flag inFlag = mFlags.registerRequired(File.class, "FILE", "SAM/BAM format files containing coordinate-sorted reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('o', CommonFlags.OUTPUT_FLAG, File.class, "FILE", "name for output SAM/BAM file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(LEGACY_CIGARS, "if set, use legacy cigars in output").setCategory(UTILITY);
    mFlags.registerOptional(X_ALTERNATE_SAM_HEADER, File.class, "FILE", "treat all SAM records as having the supplied header").setCategory(UTILITY);
    SamFilterOptions.registerMaskFlags(mFlags);
    SamFilterOptions.registerMinMapQFlag(mFlags);
    SamFilterOptions.registerMaxHitsFlag(mFlags, 'c');
    SamFilterOptions.registerMaxASMatedFlag(mFlags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(mFlags, 'u');
    SamFilterOptions.registerExcludeMatedFlag(mFlags);
    SamFilterOptions.registerExcludeUnmatedFlag(mFlags);
    SamFilterOptions.registerExcludeUnmappedFlag(mFlags);
    SamFilterOptions.registerExcludeUnplacedFlag(mFlags);
    SamFilterOptions.registerExcludeDuplicatesFlag(mFlags);
    SamFilterOptions.registerRestrictionFlag(mFlags);
    SamFilterOptions.registerBedRestrictionFlag(mFlags);
    CommonFlags.initThreadsFlag(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    CommonFlags.initNoGzip(mFlags);

    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);

    mFlags.setValidator(new SamFileMergeValidator());
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final boolean createIndex = !mFlags.isSet(CommonFlags.NO_INDEX);
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final boolean legacy = mFlags.isSet(LEGACY_CIGARS);
    final int numberThreads = CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
    final SamFilterParams filterParams = SamFilterOptions.makeFilterParamsBuilder(mFlags).create();
    File output = (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
    if ((output != null) && CommonFlags.isStdio(output)) { // Allow "-" as an alternative for writing to stdout.
      output = null;
    }
    SequencesReader template = (!mFlags.isSet(CommonFlags.TEMPLATE_FLAG)) ? null : SequencesReaderFactory.createDefaultSequencesReader((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG));

    final Collection<File> inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY).getFileList(mFlags);
    final SamMerger merger = new SamMerger(createIndex, gzip, legacy, numberThreads, filterParams, true, false);
    final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, true);
    final SAMFileHeader uberHeader;
    if (mFlags.isSet(X_ALTERNATE_SAM_HEADER)) {
      final File altHeaderFile = (File) mFlags.getValue(X_ALTERNATE_SAM_HEADER);
      uberHeader = SamUtils.getSingleHeader(altHeaderFile);
    } else {
      uberHeader = SamUtils.getUberHeader(template, inputs.getSamFiles());
    }
    merger.mergeSamFiles(inputs.getSamFiles(), inputs.getCalibrationFiles(), output, out, template, uberHeader, true, true);
    return 0;
  }

}
