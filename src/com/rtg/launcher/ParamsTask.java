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

import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Common code for those modules/tasks that can be configured via a ModuleParams.
 * @param <P> type of underlying implementation.
 * @param <S> type of statistics object.
 */
public abstract class ParamsTask<P extends ModuleParams, S extends Statistics> implements IORunnable {

  protected P mParams;  //see bug 1447 for why this is not final
  protected OutputStream mReportStream;
  protected final UsageMetric mUsageMetric;

  /** Declared here to keep the testing boffins happy */
  protected S mStatistics;

  /**
   * @param params parameters for the build and search.
   * @param reportStream stream to write statistics to.
   * @param stats statistics object to populate.
   * @param usageMetric keeps track of usage for the current module.
   */
  protected ParamsTask(final P params, final OutputStream reportStream, S stats, final UsageMetric usageMetric) {
    assert params != null;
    mParams = params;
    mReportStream = reportStream;
    mStatistics = stats;
    mUsageMetric = usageMetric;
  }

  /**
   * @return usage counter.
   */
  public long usage() {
    return mUsageMetric.getMetric();
  }

  /**
   * Execute the current task.  Does not return until the task
   * (including possible subtasks) have completed. This implementation
   * manages creation of the output directory, calling the exec
   * method, creating a done file, and closing the params upon
   * completion.
   * @throws IOException if there is a problem.
   */
  @Override
  public void run() throws IOException {
    final File dir = mParams.directory();
    //System.err.println("Making directory:" + dir.getAbsolutePath());
    if (!dir.exists() && !dir.mkdirs()) {
      throw new NoTalkbackSlimException(ErrorType.DIRECTORY_NOT_CREATED, dir.getPath());
    }
    try {
      exec();
      mStatistics.printStatistics(mReportStream);
      mStatistics.generateReport();
    } finally {
      mParams.close();
    }
  }

  /**
   * Subclasses should do all their work here
   * @throws IOException if there is a problem.
   */
  protected abstract void exec() throws IOException;

  /**
   * @return the parameters for this task.
   */
  public P parameters() {
    return mParams;
  }

  @Override
  public String toString() {
    return mParams.toString();
  }

  /**
   * Get the map containing the statistics from the run
   * @return the map
   */
  public S getStatistics() {
    return mStatistics;
  }

}
