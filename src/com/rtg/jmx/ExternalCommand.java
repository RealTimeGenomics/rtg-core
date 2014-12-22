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
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Run a command and return the lines of output. The current
 * implementation assumes that the command produces very little output
 * (such that OS buffers are sufficient to capture stdout)
 */
public class ExternalCommand {

  private static final String LS = System.lineSeparator();

  private final String[] mCommand;

  /**
   * Creates a new <code>ExternalCommand</code> instance.
   *
   * @param command the command arguments.
   */
  public ExternalCommand(String... command) {
    mCommand = command;
  }

  /**
   * Run the configured command and returns the lines of output
   *
   * @return an array of the lines of output, or null if there was a problem executing the command.
   * @exception IOException if an error occurs.
   */
  public String[] runCommand() throws IOException {
    return runCommand(mCommand);
  }


  protected String[] runCommand(String... command) throws IOException {
    try {
      final ProcessBuilder pb = new ProcessBuilder(command);
      final Process p = pb.start();
      // The command will write out less than 1kb of data, so don't need a thread to pull the output.
      p.waitFor();
      final String[] result = readAll(p.getInputStream()).split(LS);
      p.getErrorStream().close();
      p.getOutputStream().close();
      if (p.exitValue() != 0) {
        return null;
      }
      return result;
    } catch (InterruptedException ie) {
      return null;
    }
  }

  static String readAll(final InputStream is) throws IOException {
    try {
      final InputStreamReader isr = new InputStreamReader(is);

      final StringBuilder sb = new StringBuilder();
      final char[] buffer = new char[4096];
      int length;
      while ((length = isr.read(buffer)) != -1) {
        sb.append(buffer, 0, length);
      }
      return sb.toString();
    } finally {
      is.close();
    }
  }

}
