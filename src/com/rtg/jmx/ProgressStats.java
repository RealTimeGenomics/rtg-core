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
package com.rtg.jmx;

import java.io.IOException;

import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Current progress status
 */
public class ProgressStats implements MonStats {

  private static final String LS = System.lineSeparator();

  @Override
  public void addHeader(Appendable out) throws IOException {
    final String cl = CommandLine.getCommandLine();
    if (cl != null) {
      out.append("# Command line = ").append(cl).append(LS);
    }
  }

  @Override
  public void addColumnLabelsTop(Appendable out) throws IOException {
    out.append("Current ");
  }

  @Override
  public void addColumnLabelsBottom(Appendable out) throws IOException {
    out.append("progress");
  }

  @Override
  public void addColumnData(Appendable out) throws IOException {
    out.append(Diagnostic.lastProgress());
  }
}
