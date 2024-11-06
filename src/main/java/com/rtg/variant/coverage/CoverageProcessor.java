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

package com.rtg.variant.coverage;

import java.io.Closeable;
import java.io.IOException;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.util.Environment;

/**
 * Adaptor class that produces no output.
 */
@JumbleIgnore
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
  public void finalCoverageRegion(String name, int start, int end, double coverage) throws IOException { }

}
