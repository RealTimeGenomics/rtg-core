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
package com.rtg.launcher;


import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.rtg.jmx.LocalStats;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.diagnostic.Spy;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogFile;
import com.rtg.util.io.LogStream;

/**
 * Basic handling of logged command line modules. This class assumes
 * the modules have the concept of an output directory into which the
 * log is written.
 *
 */
public abstract class LoggedCli extends AbstractCli {

  static final String LOG_EXT = ".log";

  private File mDirectory = null;

  private boolean mCleanDirectory = false;


  /**
   * Determine output directory from args
   * @return output directory
   */
  protected abstract File outputDirectory();

  protected void createDirectory(File dir) {
    if (!dir.exists()) {
      mDirectory = dir;
    }
    if (!dir.isDirectory() && !dir.mkdirs()) {
      throw new NoTalkbackSlimException(ErrorType.DIRECTORY_NOT_CREATED, dir.getPath());
    }
  }

  protected void cleanDirectory() {
    if (mDirectory != null) {
      mCleanDirectory = true;
    }
  }

  protected List<DiagnosticListener> initializeOtherListeners() {
    return new ArrayList<>();
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    mDirectory = null;         // Resets state between multiple invocations
    mCleanDirectory = false;   // Resets state between multiple invocations
    final File outputDir = outputDirectory();
    createDirectory(outputDir);
    if (LocalStats.MON_DEST_OUTDIR.equals(System.getProperty(LocalStats.MON_DEST))) {
      System.setProperty(LocalStats.MON_DEST, new File(outputDir, "jmxmon.log").toString());
      LocalStats.startRecording();
    }
    final String prefix = moduleName().length() > 0 ? moduleName() : applicationName();
    final String logFile = prefix + LOG_EXT;
    final LogStream outputLog = new LogFile(new File(outputDir, logFile));
    final long startTime = System.currentTimeMillis();
    boolean successful = false;
    try {
      initializeLogs(outputLog);
      try {
        final List<DiagnosticListener> listeners = initializeOtherListeners();
        try {
          final int ret = mainExec(out, outputLog);
          if (ret == 0) {
            successful = true;
            final String time = "Finished successfully in " + timeDifference(System.currentTimeMillis(), startTime) + " s.";
            final File doneFile = new File(outputDir, "done");
            try (PrintStream done = new PrintStream(new FileOutputStream(doneFile))) {
              done.println(time);
            }

          }
          return ret;
        } finally {
          for (final DiagnosticListener list : listeners) {
            Diagnostic.removeListener(list);
          }
        }
      } catch (final SlimException e) {
        e.logException();
        throw e;
      } catch (final IOException | RuntimeException | Error e) {
        Diagnostic.userLog(e);
        throw e;
      }
    } finally {
      Spy.report();
      final String time = (successful ? "Finished successfully" : "Run failed") + " in " + timeDifference(System.currentTimeMillis(), startTime) + " s.";
      Diagnostic.userLog(time);
      Diagnostic.progress(time);
      Diagnostic.closeLog();
      if (mCleanDirectory) {
        FileUtils.deleteFiles(mDirectory);
      }
    }
  }

  /**
   * Subclasses implement this method to do the work of the
   * module. Modules can assume that logging has been set up, flags
   * are configured, and that exceptions will be handled by the
   * caller.
   *
   * @param out the <code>OutputStream</code> to be used for standard output.
   * @param log a <code>LogStream</code> to be used for logging.
   * @return main return code. 0 for usual operation, non-zero in the case of an error.
   * @exception IOException if an error occurs.
   */
  protected abstract int mainExec(OutputStream out, LogStream log) throws IOException;

  /**
   * Do common tasks for initializing logs.
   * @param initLog where to send the inital logging.
   */
  protected void initializeLogs(final LogStream initLog) {
    if (initLog != null) {
      Diagnostic.setLogStream(initLog);
      Diagnostic.logEnvironment();
      Diagnostic.progress("Started");
    } else {
      Diagnostic.setLogStream();
    }
    Diagnostic.userLog("Command line arguments: " + mFlags.getCommandLine());
    Diagnostic.userLog("Run Id: " + CommandLine.getRunId());
  }

  protected static long timeDifference(long currentTime, long startTime) {
    return (currentTime - startTime) / 1000;
  }
}
