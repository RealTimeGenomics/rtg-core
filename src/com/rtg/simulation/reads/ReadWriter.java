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
package com.rtg.simulation.reads;

import java.io.Closeable;
import java.io.IOException;

import com.rtg.reader.SdfId;

/**
 * Writes/processes reads generated from a {@link Machine}
 */
public interface ReadWriter extends Closeable {

  /**
   * Specifies the list of template sets used during generation.
   * @param templateIds an array containing an ID for each template set
   */
  void identifyTemplateSet(SdfId... templateIds);

  /**
   * Specifies the original reference template used during mutated genome generation.
   * @param referenceId the ID of the original reference template.
   */
  void identifyOriginalReference(SdfId referenceId);

  /**
   * Write a single end read
   * @param name name of read
   * @param data sequence data for read
   * @param qual quality data for read
   * @param length length of read
   * @throws IOException whenever
   */
  void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException;

  /**
   * Write left end of paired end read
   * @param name name of read
   * @param data sequence data for read
   * @param qual quality data for read
   * @param length length of read
   * @throws IOException whenever
   */
  void writeLeftRead(String name, byte[] data, byte[] qual, int length) throws IOException;

  /**
   * Write right end of paired end read
   * @param name name of read
   * @param data sequence data for read
   * @param qual quality data for read
   * @param length length of read
   * @throws IOException whenever
   */
  void writeRightRead(String name, byte[] data, byte[] qual, int length) throws IOException;

  /**
   * Returns a count of the reads written by this ReadWriter
   * @return the count of reads written
   */
  int readsWritten();
}
