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
