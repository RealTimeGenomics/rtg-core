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

import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;

/**
 */
class AlleleStatReader implements AutoCloseable, MetaSnpReader {

  private final BufferedReader mReader;
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

  @Override
  public void close() throws IOException {
    mReader.close();
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

  @Override
  public List<String> samples() {
    return mSamples;
  }

  @Override
  public MetaSnpLine nextLine() throws IOException {
    String str;
    do {
      str = mReader.readLine();
    } while (str != null && str.startsWith("#"));

    if (str == null) {
      return null;
    } else {
      return MetaSnpLine.create(str, mLineNumber++);
    }
  }
}
