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
package com.rtg.calibrate;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedUtils;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.io.IOIterator;
import com.rtg.util.io.InputFileUtils;
import com.rtg.vcf.VcfReader;

/**
 * Create calibration for quality in given SAM files
 */
public class RecalibrateCli extends AbstractCli {

  private static final String FORCE_FLAG = "force";
  private static final String COVARIATE_FLAG = "Xcovariate";
  private static final String MERGE_FLAG = "merge";
  private static final String EXCLUDE_VCF_FLAG = "exclude-vcf";
  private static final String EXCLUDE_BED_FLAG = "exclude-bed";

  @Override
  public String moduleName() {
    return "calibrate";
  }

  @Override
  public String description() {
    return "create calibration data from SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Creates quality calibration files for all supplied SAM/BAM files.");
    CommonFlagCategories.setCategories(mFlags);

    mFlags.setValidator(flags -> {
      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }
      if (!checkBedFileFlag(flags)) {
        return false;
      }
      if (!CommonFlags.validateTemplate(flags)) {
        return false;
      }

      final boolean force = flags.isSet(FORCE_FLAG);
      if (!force) {
        try {
          final List<File> files = CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, false);
          for (File f : files) {
            final File calibration = new File(f.getPath() + Recalibrate.EXTENSION);
            if (calibration.exists()) {
              flags.setParseMessage("Calibration file already exists: " + calibration.getPath());
              return false;
            }
          }
        } catch (IOException ioe) {
          Diagnostic.error("Exception reading file list: " + ioe.getMessage());
          return false;
        }
      }

      return true;
    });
    final Flag inFlag = mFlags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(CommonFlagCategories.INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    CommonFlags.initNoMaxFile(mFlags);
    CommonFlags.initThreadsFlag(mFlags);
    final Flag listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, "SDF", "SDF containing reference genome against which reads were aligned").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('f', FORCE_FLAG, "force overwriting of calibration files").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional('m', MERGE_FLAG, File.class, "file", "merge records and calibration files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    bedFileFlag(mFlags);
    // These next two will probably get included into map
    mFlags.registerOptional(EXCLUDE_VCF_FLAG, File.class, "FILE", "VCF containing sites of known variants to exclude from calibration").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(EXCLUDE_BED_FLAG, File.class, "FILE", "BED containing regions to exclude from calibration").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    final Flag covariates = mFlags.registerOptional('c', COVARIATE_FLAG, CovariateEnum.class, "COVARIATE", "covariates to recalibrate on").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    covariates.setMaxCount(Integer.MAX_VALUE);
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);
  }

  /**
   * Initialise the calibration flags for bed file and known sites
   * @param flags flags object to add bed file to
   */
  public static void bedFileFlag(CFlags flags) {
    flags.registerOptional(CommonFlags.BED_REGIONS_FLAG, File.class, "FILE", "restrict calibration to mappings falling within the supplied BED regions").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Check the bed file flag for file existence.
   * @param flags the flags to check.
   * @return true if there is nothing wrong, false otherwise.
   */
  public static boolean checkBedFileFlag(CFlags flags) {
    if (flags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
      final File bedFile = (File) flags.getValue(CommonFlags.BED_REGIONS_FLAG);
      if (!bedFile.exists()) {
        flags.setParseMessage("The specified BED file, \"" + bedFile.getPath() + "\", does not exist.");
        return false;
      } else if (bedFile.isDirectory()) {
        flags.setParseMessage("The specified BED file, \"" + bedFile.getPath() + "\", is a directory.");
        return false;
      }
    }
    return true;
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final List<File> fileCollection = InputFileUtils.removeRedundantPaths(CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false));
    final List<CovariateEnum> cs;
    if (mFlags.isSet(COVARIATE_FLAG)) {
      final List<Object> vals = mFlags.getValues(COVARIATE_FLAG);
      cs = vals.stream().map(o -> (CovariateEnum) o).collect(Collectors.toList());
    } else {
      cs = CovariateEnum.DEFAULT_COVARIATES;
    }
    final File templateSdf = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    final SequencesReader template = SequencesReaderFactory.createDefaultSequencesReader(templateSdf);
    SdfUtils.validateNoDuplicates(template, false);
    ReferenceRegions regions = mFlags.isSet(CommonFlags.BED_REGIONS_FLAG) ? BedUtils.regions((File) mFlags.getValue(CommonFlags.BED_REGIONS_FLAG)) : null;
    if (regions == null && (mFlags.isSet(EXCLUDE_VCF_FLAG) || mFlags.isSet(EXCLUDE_BED_FLAG))) {
      Diagnostic.userLog("Getting reference non-N regions");
      regions = Calibrator.getNonNRegions(template, null);
  }
    if (mFlags.isSet(EXCLUDE_VCF_FLAG)) {
      final File vcf = (File) mFlags.getValue(EXCLUDE_VCF_FLAG);
      Diagnostic.userLog("Loading exclusion VCF: " + vcf);
      try (IOIterator<? extends SequenceNameLocus> r = VcfReader.openVcfReader(vcf)) {
        regions.subtract(r);
      }
    }
    if (mFlags.isSet(EXCLUDE_BED_FLAG)) {
      final File bed = (File) mFlags.getValue(EXCLUDE_BED_FLAG);
      Diagnostic.userLog("Loading exclusion BED: " + bed);
      try (IOIterator<? extends SequenceNameLocus> r = BedReader.openBedReader(null, bed, 0)) {
        regions.subtract(r);
      }
    }
    Diagnostic.userLog("Starting calibration");
    if (regions != null) {
      ReaderUtils.validateRegions(template, regions);
    }
    try (Recalibrate r = new Recalibrate(template, regions)) {
      final File mergeFile = (File) mFlags.getValue(MERGE_FLAG);
      if (mergeFile == null) {
        r.doRecalibrate(fileCollection, cs, mFlags.isSet(FORCE_FLAG));
      } else {
        final int threads = CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
        r.doMergeRecalibrate(mergeFile, fileCollection, cs, threads, mFlags.isSet(FORCE_FLAG), true);
      }
    }
    Diagnostic.userLog("Finished calibration");
    return 0;
  }

}
