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

package com.rtg.util.io;

import java.io.FilterWriter;
import java.io.IOException;
import java.io.Writer;

import com.rtg.util.StringUtils;

/**
 * Adds methods that add a platform-specific line separator after the
 * write call, like <code>PrintStream.println()</code> but with proper
 * exceptions and no automatic flushing.
 */
public class LineWriter extends FilterWriter {

  /**
   * @param out the writer to wrap around
   */
  public LineWriter(Writer out) {
    super(out);
  }

  /**
   * Output a platform specific new line.
   * @throws IOException if the underlying writer exception
   */
  public void newLine() throws IOException {
    super.write(StringUtils.LS);
  }

  /**
   * Like the matching write method but adds a new line.
   * @param c the char
   * @throws IOException if the underlying writer exception
   */
  public void writeln(int c) throws IOException {
    super.write(c);
    newLine();
  }

  /**
   * Like the matching write method but adds a new line.
   * @param str the string
   * @throws IOException if the underlying writer exception
   */
  public void writeln(String str) throws IOException {
    super.write(str);
    newLine();
  }

  /**
   * Like the matching write method but adds a new line.
   * @param cbuf the buffer
   * @throws IOException if the underlying writer exception
   */
  public void writeln(char[] cbuf) throws IOException {
    super.write(cbuf);
    newLine();
  }

  /**
   * Like the matching write method but adds a new line.
   * @param cbuf the buffer
   * @param off offset
   * @param len length
   * @throws IOException if the underlying writer exception
   */
  public void writeln(char[] cbuf, int off, int len) throws IOException {
    super.write(cbuf, off, len);
    newLine();
  }

  /**
   * Like the matching write method but adds a new line.
   * @param str the string
   * @param off offset
   * @param len length
   * @throws IOException if the underlying writer exception
   */
  public void writeln(String str, int off, int len) throws IOException {
    super.write(str, off, len);
    newLine();
  }

  /**
   * Write an empty new line.
   * @throws IOException if the underlying writer exception
   */
  public void writeln() throws IOException {
    newLine();
  }

}
