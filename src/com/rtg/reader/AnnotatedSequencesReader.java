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
package com.rtg.reader;

/**
 * Interface to hold some methods that have no business being on the base SequencesReader interface
 */
public interface AnnotatedSequencesReader extends SequencesReader {
  /**
   * Get the SAM read group header associated with these sequences
   * Usually set at format time when formatting from SAM/BAM
   * @return the SAM read group header
   */
  String samReadGroup();

  /**
   * @return The comment stored in the SDF
   */
  String comment();

  /**
   * @return The command line parameters used to generate the SDF
   */
  String commandLine();
}
