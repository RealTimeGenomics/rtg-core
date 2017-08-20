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

import static com.rtg.vcf.VcfUtils.INFO_SVTYPE;

import java.util.ArrayList;

import com.rtg.vcf.VcfRecord;

/** Simple classification of copy number alteration types */
public enum CnaType {

  /** Deletion / loss */
  DEL,
  /** Duplication / gain */
  DUP,
  /** Not a copy number alteration */
  NONE;

  /**
   * Determines the type of copy number alteration of a VCF record
   * @param rec the record
   * @return the determined copy number alteration type
   */
  public static CnaType valueOf(final VcfRecord rec) {
    final ArrayList<String> svTypes = rec.getInfo().get(INFO_SVTYPE);
    if (svTypes == null || svTypes.size() != 1) {
      return CnaType.NONE;
    }
    switch (svTypes.get(0)) {
      case "DUP":
        return CnaType.DUP;
      case "DEL":
        return CnaType.DEL;
      default:
        return CnaType.NONE;
    }
  }
}
