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
