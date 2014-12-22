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
package com.rtg.variant;

import java.util.TreeSet;

import com.rtg.sam.ReaderRecord;

/**
 */
public interface SamToMatch {

  /**
   * Process a single alignment record by passing along to appropriate matchers (which
   * are likely to be supplied in constructors).
   * @param templateBytes the nucleotides in the template.
   * @param var alignment record to be processed.
   * @param varComplexContext record information for complex region logic in
   *          here.
   * @return true iff the record is processed without error.
   */
  boolean process(byte[] templateBytes, VariantAlignmentRecord var, TreeSet<VariantAlignmentRecord> varComplexContext);

  /**
   * Get the effective start position of the alignment record. That is, the earliest
   * reference nucleotide position that might be affected by the record. This
   * may differ from the start position specified in the alignment record in places
   * like all paths processing.
   * @param var record whose start position is to be extracted.
   * @return start position (0 based).
   */
  int start(ReaderRecord<?> var);

  /**
   * Get the effective end position of the alignment record. That is, the latest
   * reference nucleotide position that might be affected by the record. This
   * may differ from the end position specified in the alignment record in places
   * like all paths processing.
   * @param var record whose end position is to be extracted.
   * @return end position (0 based).
   */
  int end(VariantAlignmentRecord var);

}
