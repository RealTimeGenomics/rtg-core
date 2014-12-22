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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Removes unwanted Info field entries from a VCF record
 */
public class VcfInfoStripper implements VcfAnnotator {

  private final boolean mRemoveAll;
  private final boolean mKeepMode;
  private final HashSet<String> mInfos;

  /**
   * Remove all info fields from header and records
   * @param removeAll false if you don't actually want to do it for some reason.
   */
  VcfInfoStripper(boolean removeAll) {
    mRemoveAll = removeAll;
    mKeepMode = false;
    mInfos = null;
  }

  /**
   * Keep or remove a selected set of info fields from header and records
   * @param infoList the list of info field ids
   * @param keep true to keep values in the list, false to remove them
   */
  VcfInfoStripper(HashSet<String> infoList, boolean keep) {
    mRemoveAll = false;
    mKeepMode = keep;
    mInfos = infoList;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    if (mRemoveAll) {
      header.getInfoLines().clear();
      return;
    } else if (mInfos == null || mInfos.size() == 0) {
      return;
    }
    final Iterator<InfoField> it = header.getInfoLines().iterator();
    while (it.hasNext()) {
      final InfoField info = it.next();
      if (mKeepMode ^ mInfos.contains(info.getId())) {
        it.remove();
      }
    }
  }

  @Override
  public void annotate(VcfRecord rec) {
    if (mRemoveAll) {
      rec.getInfo().clear();
      return;
    } else if (mInfos == null || mInfos.size() == 0) {
      return;
    }
    final Iterator<Map.Entry<String, ArrayList<String>>> it = rec.getInfo().entrySet().iterator();
    while (it.hasNext()) {
      final Map.Entry<String, ArrayList<String>> e = it.next();
      if (mKeepMode ^ mInfos.contains(e.getKey())) {
        it.remove();
      }
    }
  }
}
