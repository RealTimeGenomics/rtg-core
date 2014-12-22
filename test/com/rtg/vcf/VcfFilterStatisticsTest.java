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

import static com.rtg.util.StringUtils.LS;

import java.util.Collections;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.vcf.VcfFilterStatistics.Stat;

import junit.framework.TestCase;

/**
 */
@SuppressWarnings("fallthrough")
public class VcfFilterStatisticsTest extends TestCase {

  public void testEnum() {
    TestUtils.testEnum(Stat.class, "[SAME_AS_REF_FILTERED_COUNT, ALL_SAME_AS_REF_FILTERED_COUNT, AMBIGOUS_FILTERED_COUNT, READ_DEPTH_FILTERED_COUNT, FAILED_KEEP_COUNT, NOT_SNP_COUNT, GENOTYPE_QUALITY_POSTERIOR_FILTERED_COUNT, QUALITY_FILTERED_COUNT, ALLELE_BALANCE_FILTERED_COUNT, DENSITY_WINDOW_COUNT, EXCLUDE_BED_COUNT, INCLUDE_BED_COUNT, WRITTEN_COUNT, TOTAL_COUNT, AVR_SCORE_FILTERED_COUNT, OVERLAP_COUNT, SNP_COUNT, DENOVO_SCORE, COMBINED_READ_DEPTH_FILTERED_COUNT]");
  }

  public void test() {
    final VcfFilterStatistics stats = new VcfFilterStatistics();
    stats.setFilterTags(Collections.singletonMap("Foo", 0));
    stats.setInfoTags(Collections.singletonMap("Bar", 0));
    final MemoryPrintStream stream = new MemoryPrintStream();
    stats.printStatistics(stream.outputStream());
    final String expected = "" + LS
        + "Total records : 0" + LS
        + "Remaining records : 0" + LS
        ;
    assertEquals(expected, stream.toString());
    stream.reset();
    stats.increment(Stat.TOTAL_COUNT);
    stats.incrementInfoTag("Bar");
    stats.incrementFilterTag("Foo");
    stats.printStatistics(stream.outputStream());
    final String expected2  = "" + LS
                            + "Total records : 1" + LS
                            + "Filtered due to Foo : 1" + LS
                            + "Filtered due to Bar : 1" + LS
                            + "Remaining records : 0" + LS
        ;
    assertEquals(expected2, stream.toString());
    stream.reset();

    for (int i = 0; i < 15; i++) {
      for (VcfFilterStatistics.Stat s : Stat.values()) {
        if (s.ordinal() > i - 2) {
          stats.increment(s);
        }
      }
    }

    stats.incrementInfoTag("Bar");
    stats.incrementFilterTag("Foo");
    stats.printStatistics(stream.outputStream());
    TestUtils.containsAll(stream.toString(),
        "Total records",
        "Filtered due to Foo",
        "Filtered due to Bar",
        "Filtered due to quality",
        "Filtered due to genotype quality",
        "Filtered due to AVR score",
        "Filtered due to sample read depth",
        "Filtered due to combined read depth",
        "Filtered due to ambiguity ratio",
        "Filtered due to allele balance",
        "Filtered due to same as reference",
        "Filtered due to all samples same as reference",
        "Filtered due to not a SNP",
        "Filtered due to simple SNP",
        "Filtered due to not in keep set",
        "Filtered due to overlap with previous",
        "Filtered due to density window",
        "Filtered due to include file",
        "Filtered due to exclude file",
        "Remaining records");
  }
}
