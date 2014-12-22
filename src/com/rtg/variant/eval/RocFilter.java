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

package com.rtg.variant.eval;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Interface defining ROC filters
 */
public enum RocFilter {

  /** accepts everything **/
  ALL {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return true;
    }
  },
  /** all homozygous **/
  HOMOZYGOUS {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return VcfUtils.isHomozygous(rec, sample);
    }
  },
  /** all heterozygous **/
  HETEROZYGOUS {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return VcfUtils.isHeterozygous(rec, sample);
    }
  },
  /** all complex calls **/
  COMPLEX {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return VcfUtils.isComplexScored(rec);
    }
  },
  /** all simple (non complex) calls **/
  SIMPLE {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return !VcfUtils.isComplexScored(rec);
    }
  },
  /** homozygous complex calls **/
  HOMOZYGOUS_COMPLEX {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return VcfUtils.isComplexScored(rec) && VcfUtils.isHomozygous(rec, sample);
    }
  },
  /** homozygous simple (non complex) calls **/
  HOMOZYGOUS_SIMPLE {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return !VcfUtils.isComplexScored(rec) && VcfUtils.isHomozygous(rec, sample);
    }
  },
  /** heterozygous complex calls **/
  HETEROZYGOUS_COMPLEX {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return VcfUtils.isComplexScored(rec) && VcfUtils.isHeterozygous(rec, sample);
    }
  },
  /** heterozygous simple (non complex) calls **/
  HETEROZYGOUS_SIMPLE {
    @Override
    boolean accept(VcfRecord rec, int sample) {
      return !VcfUtils.isComplexScored(rec) && VcfUtils.isHeterozygous(rec, sample);
    }
  };

  /**
   * Tests specified record
   * @param rec record to be tested
   * @param sample sample number
   * @return if accepted returns true, false otherwise
   */
  abstract boolean accept(VcfRecord rec, int sample);

}
