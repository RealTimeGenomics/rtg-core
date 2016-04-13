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
package com.rtg.blacklist;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.Statistics;

/**
 * Does nothing
 */
public class NullStatistics implements Statistics {
  @Override
  public void printStatistics(OutputStream reportStream) throws IOException {
  }

  @Override
  public void generateReport() throws IOException {
  }
}
