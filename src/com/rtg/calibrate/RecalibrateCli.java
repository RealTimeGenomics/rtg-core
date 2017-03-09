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

import static com.rtg.launcher.CommonFlags.BED_REGIONS_FLAG;
import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_LIST_FLAG;

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
    final Flag<File> inFlag = mFlags.registerRequired(File.class, FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(CommonFlagCategories.INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    CommonFlags.initNoMaxFile(mFlags);
    CommonFlags.initThreadsFlag(mFlags);
    final Flag<File> listFlag = mFlags.registerOptional('I', INPUT_LIST_FLAG, File.class, FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    CommonFlags.initReferenceTemplate(mFlags, true);
    mFlags.registerOptional('f', FORCE_FLAG, "force overwriting of calibration files").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional('m', MERGE_FLAG, File.class, FILE, "if set, merge records and calibration files to this output file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    bedFileFlag(mFlags);
    // These next two will probably get included into map
    mFlags.registerOptional(EXCLUDE_VCF_FLAG, File.class, FILE, "VCF containing sites of known variants to exclude from calibration").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(EXCLUDE_BED_FLAG, File.class, FILE, "BED containing regions to exclude from calibration").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    final Flag<CovariateEnum> covariates = mFlags.registerOptional('c', COVARIATE_FLAG, CovariateEnum.class, "COVARIATE", "covariates to recalibrate on").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    covariates.setMaxCount(Integer.MAX_VALUE).enableCsv();
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);

    mFlags.setValidator(flags -> CommonFlags.checkFileList(flags, INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
      && CommonFlags.validateInputFile(flags, EXCLUDE_BED_FLAG, EXCLUDE_VCF_FLAG, BED_REGIONS_FLAG)
      && (!flags.isSet(MERGE_FLAG) || flags.isSet(FORCE_FLAG) || CommonFlags.validateOutputFile(flags, (File) flags.getValue(MERGE_FLAG)))
      && CommonFlags.validateTemplate(flags));
  }

  /**
   * Initialise the calibration flags for bed file and known sites
   * @param flags flags object to add bed file to
   */
  public static void bedFileFlag(CFlags flags) {
    flags.registerOptional(BED_REGIONS_FLAG, File.class, FILE, "restrict calibration to mappings falling within the supplied BED regions").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final List<File> fileCollection = InputFileUtils.removeRedundantPaths(CommonFlags.getFileList(mFlags, INPUT_LIST_FLAG, null, false));
    final List<CovariateEnum> cs;
    if (mFlags.isSet(COVARIATE_FLAG)) {
      final List<?> vals = mFlags.getValues(COVARIATE_FLAG);
      cs = vals.stream().map(o -> (CovariateEnum) o).collect(Collectors.toList());
    } else {
      cs = CovariateEnum.DEFAULT_COVARIATES;
    }
    final File templateSdf = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    final SequencesReader template = SequencesReaderFactory.createDefaultSequencesReader(templateSdf);
    SdfUtils.validateNoDuplicates(template, false);
    ReferenceRegions regions = mFlags.isSet(BED_REGIONS_FLAG) ? BedUtils.regions((File) mFlags.getValue(BED_REGIONS_FLAG)) : null;
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
      final int threads = CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
      if (mergeFile == null) {
        r.doRecalibrate(fileCollection, cs, threads, mFlags.isSet(FORCE_FLAG));
      } else {
        r.doMergeRecalibrate(mergeFile, fileCollection, cs, threads, mFlags.isSet(FORCE_FLAG), true);
      }
    }
    Diagnostic.userLog("Finished calibration");
    return 0;
  }

}
