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

import java.util.HashSet;
import java.util.Iterator;

import com.rtg.vcf.header.FilterField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Removes unwanted filter entries from a VCF record
 */
public class VcfFilterStripper implements VcfAnnotator {

  private final boolean mRemoveAll;
  private final boolean mKeepMode;
  private final HashSet<String> mFilters;

  /**
   * Remove all filters from header and records
   * @param removeAll false if you don't actually want to do it for some reason.
   */
  VcfFilterStripper(boolean removeAll) {
    mRemoveAll = removeAll;
    mKeepMode = false;
    mFilters = null;
  }

  /**
   * Keep or remove a selected set of filters from header and records
   * @param filterList the list of filter ids
   * @param keep true to keep values in the list, false to remove them
   */
  VcfFilterStripper(HashSet<String> filterList, boolean keep) {
    mRemoveAll = false;
    mKeepMode = keep;
    mFilters = filterList;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    if (mRemoveAll) {
      header.getFilterLines().clear();
      return;
    } else if (mFilters == null || mFilters.size() == 0) {
      return;
    }
    final Iterator<FilterField> it = header.getFilterLines().iterator();
    while (it.hasNext()) {
      final FilterField filter = it.next();
      if (mKeepMode ^ mFilters.contains(filter.getId())) {
        it.remove();
      }
    }
  }

  @Override
  public void annotate(VcfRecord rec) {
    if (mRemoveAll) {
      rec.getFilters().clear();
      return;
    } else if (mFilters == null || mFilters.size() == 0) {
      return;
    }
    final Iterator<String> it = rec.getFilters().iterator();
    while (it.hasNext()) {
      final String e = it.next();
      if (mKeepMode ^ mFilters.contains(e)) {
        it.remove();
      }
    }
  }
}
