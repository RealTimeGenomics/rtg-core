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

package com.rtg.report;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.LogStream;

/**
 */
public class ReportCli extends LoggedCli {

  /** module name */
  public static final String MODULE_NAME = "report";
  /** flag to specify module */
  public static final String MODULE_FLAG = "module";
  /** flag to specify report type */
  public static final String TYPE_FLAG = "report-type";

  static class ReportFlagValidator implements Validator {
    /**
     * Checks if flags are good
     * @param flags the flags
     * @return true if good
     */
    @Override
    public boolean isValid(final CFlags flags) {
      final List<Object> inputs = flags.getAnonymousValues(0);
      for (Object in : inputs) {
        final File input = (File) in;
        if (!input.exists()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified input directory, \"" + input.getPath() + "\", does not exist.");
          return false;
        }
        if (!input.isDirectory()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified input directory, \"" + input.getPath() + "\", is not a directory.");
          return false;
        }
      }

      if (inputs.size() > ((ReportableModule) flags.getFlag(MODULE_FLAG).getValue()).getMaxInputDirectories()) {
        Diagnostic.error(ErrorType.INFO_ERROR, "Too many input directories specified for module type: " + inputs.size());
        return false;
      }

      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }

      if (flags.getFlag(TYPE_FLAG).getValue() != ReportType.HTML) {
        Diagnostic.error(ErrorType.INFO_ERROR, "Can only produce HTML reports at the moment.");
        return false;
      }

      return true;
    }
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final ReportableModule module = (ReportableModule) mFlags.getFlag(MODULE_FLAG).getValue();
    final ReportType type = (ReportType) mFlags.getFlag(TYPE_FLAG).getValue();
    final List<Object> ins = mFlags.getAnonymousValues(0);

    final File[] inputs = new File[ins.size()];
    for (int i = 0; i < ins.size(); i++) {
      inputs[i] = (File) ins.get(i);
    }

    final Report report;
    switch (module) {
      case SPECIES:
        report = new SpeciesReport();
        break;
      case SIMILARITY:
        report = new SimilarityReport();
        break;
      case MAP:
        report = new MapReport(SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnmapped(true).excludeUnplaced(true).create());
        break;
      default:
        throw new UnsupportedOperationException("Don't know about module: " + module);
    }
    report.generateReport(type, inputs, outputDirectory());

    return 0;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * set up a flags object for this module
   * @param flags the flags to set up
   */
  public void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Creates a report for outputs from RTG modules.");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired(File.class, CommonFlags.DIR, "RTG module output directory").setCategory(INPUT_OUTPUT)
      .setMinCount(1).setMaxCount(Integer.MAX_VALUE);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, CommonFlags.DIR, "directory to write report to").setCategory(INPUT_OUTPUT);
    flags.registerRequired('m', MODULE_FLAG, ReportableModule.class, CommonFlags.STRING, "module to generate report for").setCategory(REPORTING);
    flags.registerOptional('t', TYPE_FLAG, ReportType.class, CommonFlags.STRING, "type of report to generate", ReportType.HTML).setCategory(UTILITY);
    SamFilterOptions.registerMaxASMatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerMaxASUnmatedFlag(flags, SamFilterOptions.NO_SINGLE_LETTER);
    SamFilterOptions.registerExcludeUnmappedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);
    SamFilterOptions.registerExcludeMatedFlag(flags);
    flags.setValidator(new ReportFlagValidator());
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new ReportCli().mainExit(args);
  }

}
