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
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Params;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.LogStream;

/**
 * Provides a command-line binding between a ParamsTask and a ModuleParams. Oversees the construction
 * and population of an appropriate CFlags, which is then used to create a ModuleParams configuration
 * object.  This ModuleParams is then used as configuration for a ParamsTask.
 *
 * @param <P> type of underlying implementation.
 */
public abstract class ParamsCli<P extends Params> extends LoggedCli {

  /**
   * Get a task to execute with the specified parameters.
   *
   *
   * @param params the parameters for the current execution.
   * @param out where the standard output is written (depending on parameters may not be used).
   * @return the Task.
   * @throws IOException when there is an input output exception during set up.
   */
  protected abstract IORunnable task(P params, OutputStream out) throws IOException;

  /**
   * Uses flags to construct the parameters object.
   * This includes checking the flags for valid values.
   * @return the derived parameters.
   * @throws InvalidParamsException if there are errors in the values of the command line flags.
   * @throws IOException If an I/O error occurs
   */
  protected abstract P makeParams() throws InvalidParamsException, IOException;

  /**
   * Get parameters from command line, set up logging and execute task.
   * @param out where output s to be written.
   * @param initLog where to write the initial log before the command arguments have been parsed and the true log location determined.
   * @return exit code - 0 if all ok - 1 if command line arguments failed.
   * @throws IOException If an IO error occurs
   */
  @Override
  protected int mainExec(final OutputStream out, final LogStream initLog) throws IOException {
    final P localParams;
    try {
      localParams = makeParams();
    } catch (final InvalidParamsException e) {
      e.printErrorNoLog();
      //treat this as an error in the arguments passed to the process
      mFlags.error(mFlags.getInvalidFlagMsg());
      cleanDirectory();
      return 1;
    } catch (final SlimException e) {
      cleanDirectory();
      throw e;
    } catch (final RuntimeException e) {
      mFlags.error("There was an unknown error in your arguments");
      cleanDirectory();
      throw e;
    }
    try {
      Diagnostic.developerLog(localParams.toString());
      task(localParams, out).run();
      return 0;
    } finally {
      localParams.close();
    }
  }

  protected static Collection<File> checkFiles(final Collection<Object> values) throws InvalidParamsException {
    final Set<File> set = new HashSet<>();
    int filesNotFoundCount = 0;
    for (final Object obj : values) {
      final File samFile = (File) obj;
      if (!samFile.exists()) {
        filesNotFoundCount++;
        if (filesNotFoundCount <= 5) {
          Diagnostic.error(ErrorType.FILE_NOT_FOUND, samFile.getPath());
        }
      }
      set.add(samFile);
    }
    if (filesNotFoundCount > 0) {
      throw new InvalidParamsException(ErrorType.INFO_ERROR, filesNotFoundCount + " specified files were not found");
    }
    return set;
  }
}
