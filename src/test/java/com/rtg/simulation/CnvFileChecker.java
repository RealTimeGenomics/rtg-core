/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.simulation;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintStream;
import java.io.Reader;

/**
 */
public class CnvFileChecker {

  private final LineNumberReader mLin;
  private final PrintStream mErr;

  /**
   * @param err error message output
   * @param in input CNV file
   */
  public CnvFileChecker(PrintStream err, Reader in) {
    mErr = err;
    mLin = new LineNumberReader(in);
  }

  /**
   * Checks CNV file
   * @throws IOException whenever
   */
  public void check() throws IOException {
    if (!checkVersionLine()) {
      return;
    }
    if (!checkCommandLine()) {
      return;
    }
    if (!checkHeaderLine()) {
      return;
    }
    int count = 0;
    String lastLine = null;
    while (true) {
      final String line = mLin.readLine();
      if (line == null) {
        if (count == 0) {
          mErr.println("No data lines found");
        }
        break;
      }
      if (checkDataLine(lastLine, line)) {
        lastLine = line;
      }
      ++count;
    }
  }

  private int getIntField(final String s, final String line) {
    try {
      final int pos = Integer.parseInt(s);
      if (pos < 0) {
        mErr.println("Negative value: " + pos);
        mErr.println(line);
      }
      return pos;
    } catch (NumberFormatException e) {
      mErr.println("Not a valid number: " + s);
      mErr.println(line);
      return -1;
    }
  }

  private boolean checkDataLine(String lastLine, String line) {
    final String[] lineSplit = line.split("\t");
    if (lineSplit.length != 7) {
      mErr.println("Incorrect number of fields in data line");
      mErr.println(line);
      return false;
    }
    final String sequence = lineSplit[0];
    final String start = lineSplit[1];
    final String end = lineSplit[2];
    final String label = lineSplit[3];
    final String cn = lineSplit[4];
    final String bpcn = lineSplit[5];
    final String error = lineSplit[6];
    if (sequence.length() == 0 || sequence.matches("\\s")) {
      mErr.println("Invalid sequence name: " + sequence);
      mErr.println(line);
      return false;
    }

    final int startPos = getIntField(start, line);
    if (startPos < 0) {
      return false;
    }
    final int endPos = getIntField(end, line);
    if (endPos < 0) {
      return false;
    }
    if (endPos < startPos) {
      mErr.println("Negative length region: " + start + " " + end);
      return false;
    }

    boolean errorHappened = false;

    if (!"cnv".equals(label)) {
      mErr.println("Label not \"cnv\": " + label);
      errorHappened = true;
    }

    try {
      final double err = Double.parseDouble(error);
      if (err < 0.0) {
        mErr.println("Negative error: " + err);
        errorHappened = true;
      }
    } catch (NumberFormatException e) {
      mErr.println("Error not a valid number: " + error);
      errorHappened = true;
    }

    final int bcn = getIntField(bpcn, line);
    errorHappened |= bcn < 0;
    final int c = getIntField(cn, line);
    errorHappened |= c < 0;
    if (bcn > c) {
      mErr.println("bp-cn greater than cn: " + bcn + " > " + c);
      errorHappened = true;
    }

    boolean lastErrorHappened = false;
    if (lastLine != null) {
      final String[] lastLineSplit = lastLine.split("\t");
      final String lSequence = lastLineSplit[0];
      final int lStart = Integer.parseInt(lastLineSplit[1]);
      final int lEnd = Integer.parseInt(lastLineSplit[2]);
      if (lEnd - lStart == 0) {
        mErr.println("Previous length 0");
        errorHappened = true;
        lastErrorHappened = true;
      }
      if (lSequence.equals(sequence)) {
        if (lStart >= startPos) {
          mErr.println("Start out of order: " + lStart + " >= " + start);
          errorHappened = true;
          lastErrorHappened = true;
        }
        if (lEnd != startPos) {
          mErr.println("Start does not match last end: " + lStart + " + " + lEnd + " != " + start);
         errorHappened = true;
         lastErrorHappened = true;
        }
      }
    }
    if (errorHappened) {
      if (lastErrorHappened) {
        mErr.println(lastLine);
      }
      mErr.println(line);
    }
    return true;
  }

  private boolean checkVersionLine() throws IOException {
    final String str = mLin.readLine();
    if (str == null) {
      mErr.println("Empty version line - premature end of file");
      return false;
    }
    if (!str.startsWith("#Version ")) {
      mErr.println("Version line does not start with #Version");
      mErr.println(str);
      return false;
    }
    if (str.substring("#Version ".length()).trim().length() == 0) {
      mErr.println("Version line does not start have version");
      mErr.println(str);
      return true;
    }
    return true;
  }

  private boolean checkCommandLine() throws IOException {
    final String str = mLin.readLine();
    if (str == null) {
      mErr.println("Empty command line - premature end of file");
      return false;
    }
    if (!str.startsWith("#CL\t")) {
      mErr.println("Command line does not start with #CL");
      mErr.println(str);
      return false;
    }
    if (str.substring("#CL\t".length()).trim().length() == 0) {
      mErr.println("Command line does not contain command");
      mErr.println(str);
      return true;
    }
    return true;
  }

  private boolean checkHeaderLine() throws IOException {
    final String str = mLin.readLine();
    if (str == null) {
      mErr.println("Empty header line - premature end of file");
      return false;
    }
    if (!str.equals("#Seq\tstart\tend\tlabel\tcn\tbp-cn\terror")) {
      mErr.println("Header line incorrect");
      mErr.println(str);
      return false;
    }
    return true;
  }
}
