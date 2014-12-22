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
import java.util.ArrayList;
import java.util.List;

import com.rtg.bed.BedUtils;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SdfUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.InputFileUtils;

/**
 * Create calibration for quality in given SAM files
 */
public class RecalibrateCli extends AbstractCli {

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static final String FORCE_FLAG = "force";
  private static final String COVARIATE_FLAG = "Xcovariate";
  private static final String MERGE_FLAG = "merge";
  /** flag for specifying bed files */
  public static final String BED_FILE = "bed-regions";

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Creates quality calibration files for all supplied SAM/BAM files.");
    CommonFlagCategories.setCategories(flags);

    flags.setValidator(new Validator() {
      @Override
      public boolean isValid(CFlags flags) {
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
      }
    });
    final Flag inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(CommonFlagCategories.INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    CommonFlags.initNoMaxFile(flags);
    CommonFlags.initThreadsFlag(flags);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, "SDF", "SDF containing reference genome against which reads were aligned").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('f', FORCE_FLAG, "force overwriting of calibration files").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional('m', MERGE_FLAG, File.class, "file", "merge records and calibration files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    bedFileFlag(flags);
    final Flag covariates = flags.registerOptional('c', COVARIATE_FLAG, CovariateEnum.class, "COVARIATE", "covariates to recalibrate on").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    covariates.setMaxCount(Integer.MAX_VALUE);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  /**
   * Initialise the bed file flag
   * @param flags flags object to add bed file to
   * @return the new flag
   */
  public static Flag bedFileFlag(CFlags flags) {
    return flags.registerOptional(BED_FILE, File.class, "FILE", "restrict calibration to mappings falling within the supplied BED regions").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Check the bed file flag for file existence.
   * @param flags the flags to check.
   * @return true if there is nothing wrong, false otherwise.
   */
  public static boolean checkBedFileFlag(CFlags flags) {
    if (flags.isSet(BED_FILE)) {
      final File bedFile = (File) flags.getValue(BED_FILE);
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
      cs = new ArrayList<>(vals.size());
      for (Object o : vals) {
        cs.add((CovariateEnum) o);
      }
    } else {
      cs = CovariateEnum.DEFAULT_COVARIATES;
    }
    final ReferenceRegions regions = mFlags.isSet(BED_FILE) ? BedUtils.regions((File) mFlags.getValue(BED_FILE)) : null;
    final File templateSdf = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    SdfUtils.validateHasNames(templateSdf);
    try (Recalibrate r = new Recalibrate(templateSdf, regions)) {
      final File mergeFile = (File) mFlags.getValue(MERGE_FLAG);
      if (mergeFile == null) {
        r.doRecalibrate(fileCollection, cs, mFlags.isSet(FORCE_FLAG));
      } else {
        final int threads = CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
        r.doMergeRecalibrate(mergeFile, fileCollection, cs, threads, mFlags.isSet(FORCE_FLAG), true);
      }
    }
    return 0;
  }

  @Override
  public String moduleName() {
    return "calibrate";
  }

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new RecalibrateCli().mainExit(args);
  }
}
