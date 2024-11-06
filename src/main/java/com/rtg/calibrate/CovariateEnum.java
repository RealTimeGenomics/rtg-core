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
package com.rtg.calibrate;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Enum for command line input for covariates
 */
public enum CovariateEnum {

  /** read group covariate */
  READGROUP,
  /** base quality covariate */
  BASEQUALITY,
  /** template sequence covariate */
  SEQUENCE,
  /** arm of read covariate */
  ARM,
  /** machine cycle covariate (that is, position in read as sequenced) */
  MACHINECYCLE;

  /**
   * Create array of covariates from required covariate types.
   * @param vals required covariate types
   * @param header SAM header (can be null)
   * @return covariates
   */
  public static Covariate[] getCovariates(List<CovariateEnum> vals, final SAMFileHeader header) {
    final TreeSet<CovariateEnum> s = new TreeSet<>(vals);
    final Covariate[] cs = new Covariate[s.size()];
    int i = 0;
    for (final CovariateEnum type : s) {
      cs[i++] = getCovariate(header, type, 0);
    }
    return cs;
  }

  static Covariate getCovariate(final SAMFileHeader header, final CovariateEnum type, final int length) {
    switch (type) {
      case BASEQUALITY:
        return new CovariateBaseQuality();
      case READGROUP:
        if (header == null) {
          return new CovariateReadGroup();
        } else {
          final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
          if (readGroups.size() == 1) {
            return new CovariateSingleReadGroup(readGroups.get(0).getReadGroupId());
          } else {
            return new CovariateReadGroup();
          }
        }
      case SEQUENCE:
        if (header == null) {
          return new CovariateSequence();
        } else {
          return new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
        }
      case ARM:
        return new CovariateArm();
      case MACHINECYCLE:
        return new CovariateMachineCycle(length);
      default:
        throw new IllegalArgumentException();
    }
  }

  /** Default list of covariates. */
  public static final List<CovariateEnum> DEFAULT_COVARIATES = Collections.unmodifiableList(Arrays.asList(READGROUP, SEQUENCE, BASEQUALITY));

  /**
   * Return the covariate corresponding to the given name (ignoring case).
   * @param covariateName covariate to get
   * @return covariate object or null
   */
  public static CovariateEnum getCovariate(final String covariateName) {
    for (final CovariateEnum e : values()) {
      if (e.name().equalsIgnoreCase(covariateName)) {
        return e;
      }
    }
    return null;
  }
}
