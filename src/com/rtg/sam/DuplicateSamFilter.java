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
package com.rtg.sam;

import java.util.HashSet;

import htsjdk.samtools.SAMRecord;

/**
 * A SAM filter which removes duplicate entries of reads which have not been mapped as primary.
 * This is intended to use when mapping from RTG-produced SAM/BAM output (we don't output a single primary alignment
 * for all alignments, due to lack of build-on-reference)
 * Uses quite a bit of memory (data dependent).
 */
public class DuplicateSamFilter implements SamFilter {

  /* set contains hashes of read names (with left/right arm-ness taken into account as well) */
  private final HashSet<Long> mReadNameSet = new HashSet<>();

  @Override
  public boolean acceptRecord(SAMRecord rec) {
    if (rec == null || rec.getReadName() == null) {
      return false;
    }
    if (!rec.getNotPrimaryAlignmentFlag()) {
      return true;
    }
    final long hash = internalHash(rec.getReadName(), !rec.getReadPairedFlag() || rec.getFirstOfPairFlag());
    if (mReadNameSet.contains(hash)) {
      return false;
    }
    mReadNameSet.add(hash);
    return true;
  }

  private static long internalHash(final String data, boolean firstOrSingleEnd) {
    long hash = 13L;
    final char[] charArray;
    final int length = (charArray = data.toCharArray()).length;
    for (int i = 0; i < length; ++i) {
      final char c = charArray[i];
      hash = hash * 31L + c;
    }

    //shift left to make room for a 'arm' bit
    hash <<= 1;
    if (!firstOrSingleEnd) {
      hash += 1;
    }
    return hash;
  }
}
