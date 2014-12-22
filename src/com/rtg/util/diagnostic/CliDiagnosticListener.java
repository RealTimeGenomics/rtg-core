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
package com.rtg.util.diagnostic;

import java.io.PrintStream;


/**
 * A diagnostic listener suitable for command line programs.  Error and
 * warnings are simply printed on the standard error stream.
 *
 */
public class CliDiagnosticListener implements DiagnosticListener {

  private final PrintStream mErr;
  private final PrintStream mOut;

  /**
   * Listener with warning and error output to specified stream,
   * and default info output to stdout.
   * @param err where warning and error diagnostics are to be written.
   */
  public CliDiagnosticListener(final PrintStream err) {
    this(err, System.out);
  }

  /**
   * Listener with warning and error output to specified stream,
   * and info output to other specified stream.
   * @param err where warning and error diagnostics are to be written.
   * @param out where info diagnostics are to be written.
   */
  public CliDiagnosticListener(final PrintStream err, final PrintStream out) {
    mErr = err;
    mOut = out;
  }

  @Override
  public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
    if (event instanceof WarningEvent || event instanceof ErrorEvent) {
      mErr.println(event.getMessage());
      mErr.flush();
    } else if (event instanceof InformationEvent) {
      mOut.println(event.getMessage());
      mOut.flush();
    }
  }

  @Override
  public void close() { }
}

