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

import com.rtg.vcf.VcfFilterStatistics.Stat;

import junit.framework.TestCase;

/**
 */
public class VcfInfoFilterTest extends TestCase {

  public void testMinMaxIntFilter() {
    final VcfRecord withDepth = VcfReader.vcfLineToRecord("g1\t8\t.\tAAAAA\tG\t.\tPASS\tDP=50\tGT\t1/1\t0/0");
    final VcfRecord withoutDepth = VcfReader.vcfLineToRecord("g1\t8\t.\tAAAAA\tG\t.\tPASS\t.\tGT\t1/1\t0/0");
    check(withDepth, 5, 100, true);
    check(withDepth, 50, 100, true);
    check(withDepth, 5, 50, true);
    check(withDepth, 51, 100, false);
    check(withDepth, 5, 49, false);
    check(withoutDepth, 5, 100, true);
  }

  private void check(VcfRecord rec, int min, int max, boolean expected) {
    final VcfInfoFilter filter = new VcfInfoFilter.MinMaxIntFilter(new VcfFilterStatistics(), Stat.COMBINED_READ_DEPTH_FILTERED_COUNT, min, max, VcfUtils.INFO_COMBINED_DEPTH);
    assertEquals(expected, filter.accept(rec));
  }
}
