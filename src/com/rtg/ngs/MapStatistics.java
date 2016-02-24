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
package com.rtg.ngs;

import com.rtg.launcher.Statistics;

/**
 */

public interface MapStatistics extends Statistics {

  /**
   * Resets all counts to 0.
   */
  void reset();

  /**
   * Increments the count for <code>field</code> of <code>arm</code>.
   *
   * @param field field to increment
   * @param arm arm to increment
   */
  void increment(MapStatisticsField field, Arm arm);

  /**
   * Sets the count for <code>field</code> of <code>arm</code> to the given <code>value</code>.
   *
   * @param field field to set
   * @param arm arm to set
   * @param value value to be set
   */
  void set(MapStatisticsField field, Arm arm, long value);

  /**
   * Returns the value of the <code>field</code> of <code>arm</code>.
   *
   * @param field field to retrieve from
   * @param arm arm to retrieve from
   * @return value of the field
   */
  long value(MapStatisticsField field, Arm arm);

  /**
   * Returns the total value of the <code>field</code> across all arms
   *
   * @param field field to retrieve from
   * @return value of the field
   */
  long totalValue(MapStatisticsField field);

  /**
   * Returns the value of the <code>field</code> of <code>arm</code> as a percent of the total count.
   *
   * @param field field to retrieve from
   * @param arm arm to retrieve from
   * @return value of the field as a percent
   */
  double valueAsPercent(MapStatisticsField field, Arm arm);

  /**
   * Returns the total value of the <code>field</code> across all arms as a percent of the total count.
   *
   * @param field field to retrieve from
   * @return value of the field as a percent
   */
  double totalValueAsPercent(MapStatisticsField field);

  /**
   * Merges the counts of another <code>stats</code> into this statistics.
   *
   * @param stats statistics to merge
   */
  void merge(MapStatistics stats);
}
