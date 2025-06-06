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
package com.rtg.sam;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_LIST_FLAG;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.NO_HEADER;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.TEMPLATE_FLAG;
import static com.rtg.launcher.CommonFlags.THREADS_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.stream.Collectors;

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

  private static final String LEGACY_CIGARS = "legacy-cigars";
  private static final String REMOVE_DUPLICATES = "remove-duplicates";
  private static final String X_ALTERNATE_SAM_HEADER = "Xalternate-sam-header";
  private static final String X_RENAME_WITH_RG = "Xrename-read-with-rg";

  private static class SamFileMergeValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (flags.isSet(OUTPUT_FLAG)) {
        final File o = (File) flags.getValue(OUTPUT_FLAG);
        final boolean isbam = o.getName().endsWith(SamUtils.BAM_SUFFIX);
        final File output = isbam ? o : SamUtils.getZippedSamFileName(!flags.isSet(NO_GZIP), o);
        if (!CommonFlags.validateOutputFile(flags, output)) {
          return false;
        }
      }
      if (!CommonFlags.checkFileList(flags, INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
        || !SamFilterOptions.validateFilterFlags(flags, true)
        || !CommonFlags.validateThreads(flags)) {
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "sammerge";
  }

  @Override
  public String description() {
    return "merge sorted SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Merges and filters coordinate-sorted SAM/BAM files.");
    CommonFlags.initReferenceTemplate(mFlags, TEMPLATE_FLAG, false, " to use when decoding CRAM input");
    CommonFlagCategories.setCategories(mFlags);
    final Flag<File> inFlag = mFlags.registerRequired(File.class, FILE, "SAM/BAM format files containing coordinate-sorted reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = mFlags.registerOptional('I', INPUT_LIST_FLAG, File.class, FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('o', OUTPUT_FLAG, File.class, FILE, "name for output SAM/BAM file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initForce(mFlags);
    mFlags.registerOptional(LEGACY_CIGARS, "if set, produce legacy cigars (using M rather than X or =) in output").setCategory(UTILITY);
    mFlags.registerOptional(NO_HEADER, "prevent SAM/BAM header from being written").setCategory(UTILITY);
    mFlags.registerOptional(X_RENAME_WITH_RG, "rename reads by prepending with their read group ID").setCategory(UTILITY);
    mFlags.registerOptional(X_ALTERNATE_SAM_HEADER, File.class, FILE, "treat all SAM records as having the supplied header").setCategory(UTILITY);
    mFlags.registerOptional(REMOVE_DUPLICATES, "detect and remove duplicate reads based on mapping position").setCategory(SENSITIVITY_TUNING);
    SamFilterOptions.registerInvertCriteriaFlag(mFlags);
    SamFilterOptions.registerSubsampleFlags(mFlags);
    SamFilterOptions.registerMaskFlags(mFlags);
    SamFilterOptions.registerMinReadLength(mFlags);
    SamFilterOptions.registerMinMapQFlag(mFlags);
    SamFilterOptions.registerMaxHitsFlag(mFlags, 'c');
    SamFilterOptions.registerMaxASMatedFlag(mFlags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(mFlags, 'u');
    SamFilterOptions.registerSelectReadGroup(mFlags, 'r');
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
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final boolean legacy = mFlags.isSet(LEGACY_CIGARS);
    final int numberThreads = CommonFlags.parseThreads((Integer) mFlags.getValue(THREADS_FLAG));
    final SamFilterParams filterParams = SamFilterOptions.makeFilterParamsBuilder(mFlags).findAndRemoveDuplicates(mFlags.isSet(REMOVE_DUPLICATES)).create();
    final File output = mFlags.isSet(OUTPUT_FLAG) ? (File) mFlags.getValue(OUTPUT_FLAG) : new File("-");
    final SequencesReader template = (!mFlags.isSet(TEMPLATE_FLAG)) ? null : SequencesReaderFactory.createDefaultSequencesReader((File) mFlags.getValue(TEMPLATE_FLAG));

    final Collection<File> inputFiles = new CommandLineFiles(INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY).getFileList(mFlags);
    final SamMerger merger = new SamMerger(createIndex, gzip, legacy, numberThreads, filterParams, true, false);
    merger.setRenameWithRg(mFlags.isSet(X_RENAME_WITH_RG));
    final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, true);
    final SAMFileHeader uberHeader;
    if (mFlags.isSet(X_ALTERNATE_SAM_HEADER)) {
      final File altHeaderFile = (File) mFlags.getValue(X_ALTERNATE_SAM_HEADER);
      uberHeader = SamUtils.getSingleHeader(altHeaderFile);
    } else {
      uberHeader = SamUtils.getUberHeader(template, inputs.getSamFiles());
    }
    final SAMFileHeader outHeader;
    if (mFlags.isSet(NO_HEADER)) {
      outHeader = null;
    } else {
      outHeader = uberHeader.clone();
      if (filterParams.selectReadGroups() != null) {
        outHeader.setReadGroups(outHeader.getReadGroups().stream().filter(r -> filterParams.selectReadGroups().contains(r.getId()) ^ filterParams.invertFilters()).collect(Collectors.toList()));
      }
    }
    merger.mergeSamFiles(inputs.getSamFiles(), inputs.getCalibrationFiles(), output, out, template, uberHeader, outHeader, true);
    return 0;
  }

}
