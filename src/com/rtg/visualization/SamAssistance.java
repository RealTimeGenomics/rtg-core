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

package com.rtg.visualization;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

/**
 * Uses a SAM record and unfolds it into strings that can be displayed.
 */
public interface SamAssistance {

  /**
   * Use information from SAM record (including read, cigar and GS/GC fields) and <code>templateInserts</code> to construct one or more strings that can be
   * inserted into a display.
   * @param sam SAM format record.
   * @param templateStr the template with inserts already applied (this is in screen co-ordinates and is a subset of the original with modifications added).
   * @param templateBytes the template as a byte array (the original with no modifications applied)
   * @param readStart where the read starts on the template (in screen co-ordinates).
   * @param displayDots if lines match display dot
   * @return array of strings each string is padded with spaces at start (but not end). There will usually be one string or 2 in case of CG reads.
   * @throws BadSuperCigarException on bad cigar
   */
  String[] samToReads(SAMRecord sam, String templateStr, byte[] templateBytes, int readStart, boolean displayDots) throws BadSuperCigarException;
}
