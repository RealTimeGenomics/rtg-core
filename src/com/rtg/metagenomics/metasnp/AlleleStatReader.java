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
package com.rtg.metagenomics.metasnp;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;

/**
 */
class AlleleStatReader implements AutoCloseable {

  @Override
  public void close() throws IOException {
    mReader.close();
  }
  static class Line {
    private static final int NUM_FIELDS = 7;
    String mSequence;
    int mPosition;
    byte mReference;
    int[][] mCounts = new int[4][];
    Line(String str, int lineNumber) throws IOException {
      final String[] fields = StringUtils.split(str, '\t');
      if (fields.length < NUM_FIELDS) {
        throw new IOException("Error on line: " + lineNumber + " Expected at least " + NUM_FIELDS + " columns but found " + fields.length + " ");
      }
      try {
      mSequence = fields[0];
      mPosition = Integer.parseInt(fields[1]) - 1;
      mReference = (byte) DNA.valueOf(fields[2]).ordinal();
      for (int i = 0; i < 4; i++) {
        final String count = fields[3 + i];
        final String[] samples = StringUtils.split(count, ',');
        mCounts[i] = new int[samples.length];
        for (int j = 0; j < samples.length; j++) {
          mCounts[i][j] = Integer.parseInt(samples[j]);
        }
      }
      } catch (NumberFormatException e) {
        throw new IOException("Error on line: " + lineNumber + " " + e.getMessage());
      }

    }
  }

  final BufferedReader mReader;
  private int mLineNumber = 1;
  private List<String> mSamples = null;

  AlleleStatReader(File in) throws IOException {
    this(FileUtils.createInputStream(in, false));
  }

  AlleleStatReader(InputStream in) throws IOException {
    mReader = new BufferedReader(new InputStreamReader(in));
    // First line is an uncommented header to placate R and it's column naming
    String line;
    while ((line = mReader.readLine()) != null && line.startsWith("#")) {
      if (line.startsWith("##Samples=")) {
        mSamples = parseSamples(line);
      }
    }
  }

  /**
   * Extract comma separated list of sample ids
   * @param line line from parsed file
   * @return sample ids as a list
   */
  private List<String> parseSamples(String line) {
    final String sampleString = line.substring(line.indexOf('=') + 1);
    final String[] parts = StringUtils.split(sampleString, ',');
    return Arrays.asList(parts);
  }

  /**
   * @return the list of sample names from the header or null if the header is missing
   */
  public List<String> samples() {
    return mSamples;
  }

  Line nextLine() throws IOException {
    String str;
    do {
      str = mReader.readLine();
    } while(str != null && str.startsWith("#"));

    if (str == null) {
      return null;
    } else {
      return new Line(str, mLineNumber++);
    }
  }
}
