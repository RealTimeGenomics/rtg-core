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

import java.util.Arrays;
import java.util.Locale;

import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.annotation.AbstractDerivedAnnotation;
import com.rtg.vcf.annotation.DerivedAnnotations;

/**
 * Enumeration of possible ROC score fields
 */
public enum RocScoreField {

  /** use <code>QUAL</code> VCF field to sort */
  QUAL {
    @Override
    RocSortValueExtractor getExtractor(final String fieldName, final RocSortOrder order) {
      return new RocSortValueExtractor() {
        @Override
        double getSortValue(VcfRecord rec, int sampleNo) {
          final String qualStr = rec.getQuality();
          if (qualStr.equals(VcfRecord.MISSING)) {
            return Double.NaN;
          }
          try {
            return Double.parseDouble(qualStr);
          } catch (NumberFormatException e) {
            throw new NoTalkbackSlimException("Invalid QUAL value = " + rec.getQuality());
          }
        }
      };
    }
  },
  /** use VCF info field to provide sort value */
  INFO {
    @Override
    RocSortValueExtractor getExtractor(final String fieldName, final RocSortOrder order) {
      return new RocSortValueExtractor() {
        @Override
        double getSortValue(VcfRecord rec, int sampleNo) {
          final double val = VcfUtils.getDoubleInfoFieldFromRecord(rec, fieldName);
          if (MathUtils.approxEquals(val, 0, 0.00000001)) {
            return 0;
          }
          return val;
        }
      };
    }
  },
  /** use VCF format field to provide sort value */
  FORMAT {
    @Override
    RocSortValueExtractor getExtractor(final String fieldName, final RocSortOrder order) {
      return new RocSortValueExtractor() {
        @Override
        double getSortValue(VcfRecord rec, int sampleNo) {
          final double val = VcfUtils.getDoubleFormatFieldFromRecord(rec, sampleNo, fieldName);
          if (MathUtils.approxEquals(val, 0, 0.00000001)) {
            return 0;
          }
          return val;
        }
      };
    }
  },
  /** use derived field to provide sort value */
  DERIVED {
    @Override
    RocSortValueExtractor getExtractor(final String fieldName, final RocSortOrder order) {
      final DerivedAnnotations derived;
      final String field = fieldName.toUpperCase(Locale.getDefault());
      try {
        derived = DerivedAnnotations.valueOf(field);
      } catch (IllegalArgumentException e) {
        throw new NoTalkbackSlimException("Unrecognized derived annotation \"" + field + "\", must be one of " + Arrays.toString(DerivedAnnotations.values()));
      }
      final AbstractDerivedAnnotation anno = derived.getAnnotation();
      return new RocSortValueExtractor() {
        @Override
        double getSortValue(VcfRecord rec, int sampleNo) {
          final Double val = (Double) anno.getValue(rec, sampleNo);
          if (val == null) {
            return Double.NaN;
          }
          if (MathUtils.approxEquals(val, 0, 0.00000001)) {
            return 0;
          }
          return val;
        }
      };
    }
  };

  abstract RocSortValueExtractor getExtractor(String fieldName, RocSortOrder order);

}
