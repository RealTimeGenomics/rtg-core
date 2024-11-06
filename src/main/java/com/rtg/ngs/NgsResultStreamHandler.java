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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.ResultStreamHandler;

/**
 * Handles result streams
 */
public class NgsResultStreamHandler extends ResultStreamHandler {

  static final String OUT_SUFFIX = "out";

  /**
   * Constructor
   * @param outDir output directory
   * @param zip whether to zip output
   */
  public NgsResultStreamHandler(File outDir, boolean zip) {
    super(outDir, null, zip);
  }

  /**
   * Constructor
   * @param outDir output directory
   * @param tempFilesDir optional directory for temporary files
   * @param zip whether to zip output
   */
  public NgsResultStreamHandler(File outDir, File tempFilesDir, boolean zip) {
    super(outDir, tempFilesDir, zip);
  }

  /**
   * Get a stream to the output file.
   * @return the stream.
   * @throws IOException if the stream cannot be created
   */
  public OutputStream outStream() throws IOException {
    return createFileStream(OUT_SUFFIX);
  }

  /**
   * Get a stream to the SAM alignments file.
   * @return the stream.
   * @throws IOException if the stream cannot be created
   */
  public OutputStream matedSamStream() throws IOException {
    return createFileStream(NgsOutputParams.MATED_SAM_FILE_NAME);
  }

  /**
   * Stream to write repeats to.
   * @return the stream
   * @throws IOException if the stream cannot be created
   */
  public OutputStream repeatStream() throws IOException {
    return createFileStream(NgsOutputParams.REPEATS_FILE_NAME);
  }

  /**
   * Stream to write unmapped reads to.
   * @return the stream
   * @throws IOException if the stream cannot be created
   */
  public OutputStream unmappedStream() throws IOException {
    return createFileStream(NgsOutputParams.UNMAPPED_FILE_NAME);
  }

}
