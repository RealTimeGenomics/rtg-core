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
  /** read position covariate */
  READPOSITION,
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
    final TreeSet<CovariateEnum> s = new TreeSet<>();
    for (final CovariateEnum e : vals) {
      s.add(e);
    }
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
      case READPOSITION:
        return new CovariateReadPos(length);
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
