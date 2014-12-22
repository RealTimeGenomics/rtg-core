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
package com.rtg.graph;

/**
 * Interface to separate loading code from UI progress code
 */
public interface ProgressDelegate {

  /**
   * Will be called with number of lines processed so far, every 100 lines
   * @param progress number of lines read so far
   */
  void setProgress(int progress);

  /**
   * Will be called once per file with the number of lines processed in that file
   * @param numberLines number of lines in file
   */
  void addFile(int numberLines);

  /**
   * Will be called once when all files are loaded
   */
  void done();
}
