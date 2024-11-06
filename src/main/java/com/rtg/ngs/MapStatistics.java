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
package com.rtg.ngs;

import com.rtg.launcher.Statistics;
import com.rtg.reader.Arm;

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

}
