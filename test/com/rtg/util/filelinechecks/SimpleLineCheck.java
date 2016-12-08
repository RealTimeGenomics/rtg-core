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
package com.rtg.util.filelinechecks;

import java.util.Arrays;

import com.rtg.util.StringUtils;
import com.rtg.util.integrity.Exam;

/**
 */
public class SimpleLineCheck {

  /** For checking tab separated lines. */
  public static final SimpleLineCheck TAB_CHECK = new SimpleLineCheck("\\t");

  private final String mSeparator;

  /**
   * @param separator regular expression used to separate fields.
   */
  public SimpleLineCheck(final String separator) {
    mSeparator = separator;
  }

  /**
   * Check each field in a line.
   * @param line0 to be checked.
   * @param checks one for each field.
   * @return true - so that can be called from assert statements.
   * @throws RuntimeException if there is an error.
   */
  public boolean lineCheck(final String line0, final Object... checks) {
    //strip any trailing LS
    final String line = line0.replace(StringUtils.LS, "");
    if (line == null) {
      throw new RuntimeException("Line is null.");
    }
    final String[] split = line.split(mSeparator);
    if (split.length > checks.length) {
      throw new RuntimeException("Incorrect number of fields. line=" + StringUtils.display(line) + " checks=" + Arrays.toString(checks) + " split=" + Arrays.toString(split));
    }
    try {
      for (int i = 0; i < checks.length; ++i) {
        final Object obj = checks[i];
        if (i >= split.length && obj == null) {
          continue;
        }
        final String field = split[i];
        checkField(field, obj);
      }
    } catch (final RuntimeException e) {
      throw new RuntimeException("Error in line=" + line + " checks=" + Arrays.toString(checks) + " split=" + Arrays.toString(split), e);
    }
    return true;
  }

  private void checkField(final String field, final Object check) {
    if (check == null) {
      return;
    }
    if (check instanceof String) {
      Exam.assertEquals(check, field);
    }
    if (check instanceof Check) {
      final Check ch = (Check) check;
      ch.check(field);
      return;
    }
    if (check instanceof Double) {
      final double d1 = Double.valueOf(field);
      final double d2 = (Double) check;
      if (Double.isInfinite(d1) || Double.isInfinite(d2)) {
        Exam.assertTrue(d1 > 0 == d2 > 0); //do it this tricky way because findbugs objects to a direct equality check
      } else {
        Exam.assertEquals(d1, d2, 0.001);
      }
    }
  }

}
