/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.coverage;

import java.io.Closeable;
import java.io.IOException;

import com.rtg.util.Environment;

/**
 * Adaptor class that produces no output.
 */
public abstract class CoverageProcessor implements Closeable {

  /** Coverage file version string */
  public static final String VERSION_STRING = "#Version " + Environment.getVersion();

  static final String COVERAGE_OUTPUT_VERSION = "v1.1";

  protected static final String TB = "\t";

  /**
   * Performs any required initialisation
   * @throws IOException because it might do some output
   */
  public void init() throws IOException { }

  /**
   * @param name sequence name
   * @param position position on sequence
   * @param ih1 count of <code>IH=1</code> records
   * @param ihgt1 count of <code>IH&gt;1</code> records
   * @param coverage coverage at position
   * @throws IOException if an IO error occurs
   */
  public void finalCoveragePosition(String name, int position, int ih1, int ihgt1, double coverage) throws IOException { }

  /**
   * Sets the label to be used for the next region
   * @param label the label
   */
  public void setRegionLabel(String label) { }

  /**
   *
   * @param name sequence name
   * @param start start position
   * @param end end position
   * @param coverage coverage for the region
   * @throws IOException because it might do some output
   */
  public void finalCoverageRegion(String name, int start, int end, int coverage) throws IOException { }

}
