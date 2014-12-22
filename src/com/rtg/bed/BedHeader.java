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

package com.rtg.bed;

import java.util.ArrayList;
import java.util.List;

/**
 * BED header container.
 */
public class BedHeader {

  private final List<String> mLines = new ArrayList<>();

  /**
   * Add a BED header line.
   * @param line BED header line.
   */
  public void addHeaderLine(String line) {
    mLines.add(line);
  }

  /**
   * Get the BED header lines.
   * @return the BED header lines.
   */
  public String[] getHeaderLines() {
    return mLines.toArray(new String[mLines.size()]);
  }
}
