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
 * Specifies the interface for retrieving sequences specific for Complete Genomics data.
 *
 *
 */
public interface CgSequenceDataSource extends SequenceDataSource {

  /**
   * return the largest number of Ns seen in any of the arms
   * @return maximum N count
   */
  int getMaxNCount();

  /**
   * Method to get how many reads were skipped due to filter conditions.
   * @return number of skipped reads
   */
  int getSkippedReads();

  /**
   * Method to get the number of residues skipped due to filter conditions.
   * @return number of residues skipped
   */
  long getSkippedResidues();

  /**
   * @return the minimum read length
   */
  @Override
  long getMinLength();

  /**
   * @return the maximum read length
   */
  @Override
  long getMaxLength();
}

