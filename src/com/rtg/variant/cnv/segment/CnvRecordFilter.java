/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.segment;

import java.util.Collection;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Only retains valid SV DUP/DEL records on the chromosomes of interest
 */
public class CnvRecordFilter implements VcfFilter {

  private final Collection<String> mChrs;
  private final boolean mFilterOverlap;
  private String mLastSeq = null;
  private int mLastEnd = 0;

  /**
   * Constructor
   * @param chrs the set of chromosomes of interest
   * @param filterOverlap true if only the first of overlapping variants should be kept
   */
  public CnvRecordFilter(Collection<String> chrs, boolean filterOverlap) {
    mChrs = chrs;
    mFilterOverlap = filterOverlap;
  }

  @Override
  public void setHeader(VcfHeader header) { }

  @Override
  public boolean accept(VcfRecord rec) {
    final CnaType status = CnaType.valueOf(rec);
    if (status != CnaType.DEL && status != CnaType.DUP) { // Non SV or not a DEL/DUP
      return false;
    }
    final Integer end = VcfUtils.getIntegerInfoFieldFromRecord(rec, VcfUtils.INFO_END);
    if (end == null) {
      Diagnostic.warning("Skipping SV record without a defined END: " + rec);
      return false;
    }
    if (mFilterOverlap && rec.getSequenceName().equals(mLastSeq) && rec.getStart() + 1 < mLastEnd) { // +1 since spec says SV records include the base before the event
      Diagnostic.warning("Skipping SV record that overlaps a previous SV variant: " + rec);
      return false;
    }
    if (!mChrs.contains(rec.getSequenceName())) {
      return false;
    }
    mLastSeq = rec.getSequenceName();
    mLastEnd = end;
    return true;
  }

}
