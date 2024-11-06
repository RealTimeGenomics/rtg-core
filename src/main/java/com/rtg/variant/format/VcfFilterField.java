/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.format;

import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantParams;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 * Enum of VCF FILTER field implementations
 */
public enum VcfFilterField {

  /** Coverage Filter */
  OC(VariantFilter.COVERAGE) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "Coverage threshold exceeded");
    }
  },
  /** Ambiguity Filter */
  A(VariantFilter.AMBIGUITY) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(params), "Ambiguity exceeded " + threshold(params));
    }
    @Override
    protected String name(VariantParams params) {
      return "a" + threshold(params);
    }
    private double threshold(VariantParams params) {
      if (params != null) {
        return params.maxAmbiguity() * 100;
      }
      return DUMMY_AMBIGUITY_THRESHOLD;
    }
  },
  /** Complex Filter */
  RC(VariantFilter.FAILED_COMPLEX) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "RTG variant is a complex region");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params) {
      if (!call.isFiltered(VariantFilter.HYPER_COMPLEX)) {
        super.updateRecord(rec, call, params);
      }
    }
  },
  /** Hypercomplex Filter */
  RX(VariantFilter.HYPER_COMPLEX) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "RTG variant contained within hypercomplex region");
    }
  },
  /** Equivalent Variant Filter */
  RCEQUIV(VariantFilter.COMPLEX_EQUIVALENT) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "RTG variant is equivalent to the previous variant");
    }
  },
  /** IonTorrent Filter */
  IONT(VariantFilter.IONTORRENT) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "IonTorrent specific filter applied");
    }
  },
  /** BED Filter */
  BED(VariantFilter.BED_REGION) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "Variant falls outside the specified target regions");
    }
  },

  /** Pass */
  PASS(null) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) { }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params) {
      if (!call.isFiltered()) {
        rec.addFilter(name());
      }
    }
  },
  /**
   *  Other Filter
   *  Note: This filter should be called after all other normal filters have been checked.
   *  That means that this one should always be the last enum value.
   */
  OTHER(VariantFilter.OTHER) {
    @Override
    public void updateHeader(VcfHeader header, VariantParams params) {
      header.addFilterField(name(), "Variant is invalid for unknown reasons");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params) {
      if (rec.getFilters().isEmpty() && call.isFiltered()) {
        rec.addFilter(name());
      }
    }
  };

  private static final double DUMMY_AMBIGUITY_THRESHOLD = -1.0;

  private final VariantFilter mFilter;

  /**
   * @param filter the filter that is being used as a trigger.
   */
  VcfFilterField(VariantFilter filter) {
    mFilter = filter;
  }

  /**
   * Update the VCF header with the field description.
   * @param header the VCF header for which the field description will be added.
   * @param params the variant output options params.
   */
  public abstract void updateHeader(VcfHeader header, VariantParams params);

  /**
   * Update the VCF record with the value for the FILTER field.
   * @param rec the VCF record for which the FILTER field value will be added.
   * @param call the variant data to use to update the VCF record.
   * @param params the variant output options params.
   */
  public void updateRecord(VcfRecord rec, Variant call, VariantParams params) {
    if (call.isFiltered(mFilter)) {
      rec.addFilter(name(params));
    }
  }

  /**
   * Get the filter name.
   * @param params the params to use in the construction of the name.
   * @return the name of the filter.
   */
  protected String name(VariantParams params) {
    return name();
  }

}
