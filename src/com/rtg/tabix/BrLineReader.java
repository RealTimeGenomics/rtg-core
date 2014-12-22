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
package com.rtg.tabix;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * Makes a LineReader around a BufferedReader so we can deal with them interchangeably with TabixLineReader
 */
public class BrLineReader implements LineReader {

  private final BufferedReader mBr;

  /**
   * Makes a LineReader around a BufferedReader
   * @param br the inner reader
   */
  public BrLineReader(BufferedReader br) {
    mBr = br;
  }

  @Override
  public String readLine() throws IOException {
    return mBr.readLine();
  }

  @Override
  public void close() throws IOException {
    mBr.close();
  }
}
