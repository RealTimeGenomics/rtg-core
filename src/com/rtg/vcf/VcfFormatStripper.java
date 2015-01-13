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
package com.rtg.vcf;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Removes unwanted Format field entries from a VCF record
 */
public class VcfFormatStripper implements VcfAnnotator {

  private final boolean mKeepMode;
  private final Set<String> mFormats;

  private boolean mKeepRecord = true;

  /**
   * Keep or remove a selected set of format fields from header and records
   * @param formatList the list of format field ids
   * @param keep true to keep values in the list, false to remove them
   */
  VcfFormatStripper(Set<String> formatList, boolean keep) {
    mKeepMode = keep;
    mFormats = formatList;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    if (mFormats == null || mFormats.size() == 0) {
      return;
    }
    final Iterator<FormatField> it = header.getFormatLines().iterator();
    while (it.hasNext()) {
      final FormatField format = it.next();
      if (mKeepMode ^ mFormats.contains(format.getId())) {
        it.remove();
      }
    }
  }

  @Override
  public void annotate(VcfRecord rec) {
    if (mFormats == null || mFormats.size() == 0) {
      return;
    }
    final Iterator<Map.Entry<String, ArrayList<String>>> it = rec.getFormatAndSample().entrySet().iterator();
    while (it.hasNext()) {
      final Map.Entry<String, ArrayList<String>> e = it.next();
      if (mKeepMode ^ mFormats.contains(e.getKey())) {
        it.remove();
      }
    }
    if (rec.getFormatAndSample().size() == 0) {
      mKeepRecord = false;
    }
  }

  boolean keepRecord() {
    return mKeepRecord;
  }
}
