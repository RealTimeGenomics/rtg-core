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

package com.rtg.metagenomics;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.metagenomics.MetagenomicsWrapperCli.Platform;
import com.rtg.ngs.AbstractSdfOutputProcessor;
import com.rtg.ngs.MapCli;
import com.rtg.ngs.MapFCli;
import com.rtg.protein.MapXCli;
import com.rtg.reader.FormatCli;
import com.rtg.reader.ReaderUtils;
import com.rtg.report.CombinedReport;
import com.rtg.sam.SamCommandHelper;
import com.rtg.usage.UsageMetric;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogStream;

/**
 * The task code for the metagenomics wrapper.
 */
class MetagenomicsWrapperTask extends ParamsTask<MetaPipelineParams, NoStatistics> {

  private final LogStream mLogStream;
  private final PrintStream mErr;
  private final List<File> mOutputDirectories = new ArrayList<>();

  private int mReturnCode = 0;

  private static class NoListenerMapF extends MapFCli {
    @Override
    protected CliDiagnosticListener initializeMainListener(PrintStream err, PrintStream out) {
      // Prevent double printing of messages.
      return new CliDiagnosticListener(NullStreamUtils.getNullPrintStream(), NullStreamUtils.getNullPrintStream());
    }
  }

  private static class NoListenerMap extends MapCli {
    @Override
    protected CliDiagnosticListener initializeMainListener(PrintStream err, PrintStream out) {
      // Prevent double printing of messages.
      return new CliDiagnosticListener(NullStreamUtils.getNullPrintStream(), NullStreamUtils.getNullPrintStream());
    }
  }

  private static class NoListenerSpecies extends SpeciesCli {
    @Override
    protected CliDiagnosticListener initializeMainListener(PrintStream err, PrintStream out) {
      // Prevent double printing of messages.
      return new CliDiagnosticListener(NullStreamUtils.getNullPrintStream(), NullStreamUtils.getNullPrintStream());
    }
  }

  private static class NoListenerMapX extends MapXCli {
    @Override
    protected CliDiagnosticListener initializeMainListener(PrintStream err, PrintStream out) {
      // Prevent double printing of messages.
      return new CliDiagnosticListener(NullStreamUtils.getNullPrintStream(), NullStreamUtils.getNullPrintStream());
    }
  }

  protected MetagenomicsWrapperTask(MetaPipelineParams params, OutputStream reportStream, UsageMetric usageMetric, LogStream logStream, PrintStream err) {
    super(params, reportStream, new NoStatistics(), usageMetric);
    mLogStream = logStream;
    mErr = err;
  }

  protected int returnCode() {
    return mReturnCode;
  }

  @Override
  protected void exec() throws IOException {
    final File output = mParams.directory();
    final File mapfOutput = new File(output, "mapf");
    final File filter = mParams.filterSdf();
    final List<String> inputFlags = new ArrayList<>();
    if (mParams.inputFile() != null) {
      final File inputFile =  mParams.inputFile();
      if (ReaderUtils.isSDF(inputFile)) {
        inputFlags.add("--" + CommonFlags.READS_FLAG);
        inputFlags.add(inputFile.getPath());
      } else {
        inputFlags.add("--format");
        inputFlags.add("fastq");
        inputFlags.add("--quality-format");
        inputFlags.add("sanger");
        inputFlags.add("--" + CommonFlags.READS_FLAG);
        inputFlags.add(inputFile.getPath());
      }
    } else {
      inputFlags.add("--format");
      inputFlags.add("fastq");
      inputFlags.add("--quality-format");
      inputFlags.add("sanger");
      final File leftFile =  mParams.inputLeft();
      final File rightFile =  mParams.inputRight();
      inputFlags.add("--" + FormatCli.LEFT_FILE_FLAG);
      inputFlags.add(leftFile.getPath());
      inputFlags.add("--" + FormatCli.RIGHT_FILE_FLAG);
      inputFlags.add(rightFile.getPath());
    }
    final List<String> mapInput = new ArrayList<>();
    if (mParams.filterSdf() != null) {
      final List<String> mapfArgs = new ArrayList<>();
      mapfArgs.addAll(Arrays.asList(
          "--" + CommonFlags.OUTPUT_FLAG, mapfOutput.getPath()
          , "--" + CommonFlags.TEMPLATE_FLAG, filter.getPath()
          // Because mapf uses the --sam-rg flag to detect platform for internal aligner config :(
          , "--" + SamCommandHelper.SAM_RG, "@RG\\tPL:" + mParams.inputPlatform() + "\\tSM:sample\\tID:id"
      ));
      mapfArgs.addAll(inputFlags);
      final MapFCli mapF = new NoListenerMapF();
      final int mapFResult = runCommand(mReportStream, mapF, mapfArgs);
      mOutputDirectories.add(mapfOutput);
      if (mapFResult != 0) {
        mReturnCode = mapFResult;
        return;
      }
      mapInput.addAll(Arrays.asList(
          "--" + CommonFlags.READS_FLAG, new File(mapfOutput, AbstractSdfOutputProcessor.UNMAPPED_SDF_FILE).getPath()
      ));
    } else {
      mapInput.addAll(inputFlags);
    }
    if (mParams.speciesSdf() != null) {
      final int result = species(mReportStream, output, mapInput);
      if (result != 0) {
        mReturnCode = result;
        return;
      }
    }
    if (mParams.proteinSdf() != null) {
      final int result = mapX(mapfOutput, mReportStream, output);
      if (result != 0) {
        mReturnCode = result;
        return;
      }
    }
    final File reportDir = new File(output, "report");
    if (!reportDir.exists() && !reportDir.mkdir()) {
      throw new IOException("Couldn't create directory: '" + reportDir + "'");
    }
    final CombinedReport report = new CombinedReport(mOutputDirectories, reportDir);
    report.makeReport();
  }

  protected static String errorRate(Platform inputMachineType) {
    return inputMachineType == Platform.ILLUMINA ? "10%" : "15%";
  }

  int mapX(File mapfOutput, OutputStream out, File output) throws IOException {
    final List<File> mapxInput = new ArrayList<>();
    if (mParams.filterSdf() != null) {
      final File filterOutput = new File(mapfOutput, AbstractSdfOutputProcessor.UNMAPPED_SDF_FILE);
      if (ReaderUtils.isPairedEndDirectory(filterOutput)) {
        mapxInput.add(new File(filterOutput, "left"));
        mapxInput.add(new File(filterOutput, "right"));
      } else {
        mapxInput.add(filterOutput);
      }
    } else {
      if (mParams.inputFile() != null) {
        final File inputFile =  mParams.inputFile();
        if (ReaderUtils.isSDF(inputFile) && ReaderUtils.isPairedEndDirectory(inputFile)) {
          mapxInput.add(new File(inputFile, "left"));
          mapxInput.add(new File(inputFile, "right"));
        } else {
          mapxInput.add(inputFile);
        }
      } else {
        final File leftFile =  mParams.inputLeft();
        final File rightFile =  mParams.inputRight();
        mapxInput.add(leftFile);
        mapxInput.add(rightFile);
      }
    }
    for (int i = 0; i < mapxInput.size(); i++) {
      final File mapxReads = mapxInput.get(i);
      final List<String> mapxFlags = new ArrayList<>();
      final File mapxOutput = new File(output, "mapx" + (i + 1));
      mOutputDirectories.add(mapxOutput);
      mapxFlags.addAll(Arrays.asList(
          "--" + CommonFlags.OUTPUT_FLAG, mapxOutput.getPath()
          , "--" + CommonFlags.TEMPLATE_FLAG, mParams.proteinSdf().getPath()
          , "--" + MapXCli.MAX_ALIGNMENT_SCORE, errorRate(mParams.inputPlatform())
          , "--" + CommonFlags.MAX_TOP_RESULTS_FLAG, "10"
      ));
      mapxFlags.add("--" + CommonFlags.READS_FLAG);
      mapxFlags.add(mapxReads.getPath());
      if (!ReaderUtils.isSDF(mapxReads)) {
        mapxFlags.add("--" + FormatCli.FORMAT_FLAG);
        mapxFlags.add("fastq");
      }
      final MapXCli mapX = new NoListenerMapX();
      final int mapXResult = runCommand(out, mapX, mapxFlags);
      if (mapXResult != 0) {
        return mapXResult;
      }
    }
    return 0;
  }

  private static class BamFilter implements FilenameFilter {
    @Override
    public boolean accept(File dir, String name) {
      return name.endsWith(".bam");
    }
  }

  int species(OutputStream out, File output, List<String> mapInput) throws IOException {
    final File mapOutput = new File(output, "map");
    final List<String> mapArgs = new ArrayList<>();
    mapArgs.addAll(Arrays.asList(
        "--" + CommonFlags.OUTPUT_FLAG, mapOutput.getPath()
        , "--" + CommonFlags.TEMPLATE_FLAG, mParams.speciesSdf().getPath()
        , "--" + CommonFlags.MAX_ALIGNMENT_MISMATCHES, errorRate(mParams.inputPlatform())
        , "--" + CommonFlags.MAX_TOP_RESULTS_FLAG, "100"
        , "--" + SamCommandHelper.SAM_RG, "@RG\\tPL:" + mParams.inputPlatform() + "\\tSM:sample\\tID:id"
    ));
    mapArgs.addAll(mapInput);

    mOutputDirectories.add(mapOutput);
    final MapCli ramMap = new NoListenerMap();
    final int mapResult = runCommand(out, ramMap, mapArgs);
    if (mapResult != 0) {
      return mapResult;
    }
    final List<String> speciesFlags = new ArrayList<>();
    final File speciesOutput = new File(output, "species");
    mOutputDirectories.add(speciesOutput);
    speciesFlags.addAll(Arrays.asList(
        "--" + CommonFlags.OUTPUT_FLAG, speciesOutput.getPath()
        , "--" + SpeciesCli.TEMPLATE_FLAG, mParams.speciesSdf().getPath()
    ));
    final String[] bamFiles = mapOutput.list(new BamFilter());
    if (bamFiles != null) {
      Arrays.sort(bamFiles);
      for (String f : bamFiles) {
        speciesFlags.add(new File(mapOutput, f).getPath());
      }
    }
    final SpeciesCli speciesCli = new NoListenerSpecies();
    return runCommand(out, speciesCli, speciesFlags);
  }

  /**
   * Execute a module with the provided arguments
   * @param out standard output destination
   * @param module a module to run
   * @param args command line arguments
   * @return 0 on success
   * @throws IOException when everything breaks
   */
  public int runCommand(OutputStream out, LoggedCli module, List<String> args) throws IOException {
    final StringBuilder command = new StringBuilder();
    try {
      command
          .append(module.applicationName())
          .append(" ")
          .append(module.moduleName())
          .append(" ");
      for (String s : args) {
        command.append(" ");
        if (s.matches(".*\\s.*")) {
          command.append("\"")
              .append(s)
              .append("\"");
        } else {
          command.append(s);
        }

      }
      out.write(("## " + command.toString()).getBytes());
      out.write(StringUtils.LS.getBytes());
    } finally {
      Diagnostic.setLogStream(mLogStream);
    }
    final List<String> commandLine = new ArrayList<>();
    commandLine.add(module.applicationName());
    commandLine.add(module.moduleName());
    commandLine.addAll(args);
    CommandLine.setCommandArgs(commandLine.toArray(new String[commandLine.size()]));
    return module.mainInit(args.toArray(new String[args.size()]), out, mErr);
  }

}
