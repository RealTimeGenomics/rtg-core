/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.multisample.cancer;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 * @author Sean A. Irvine
 */
public class PriorBedRangeLoaderTest extends TestCase {

  public void test() throws IOException {
    final File f = FileUtils.stringToFile("chr1\t1\t100\t0.42", FileHelper.createTempFile());
    try {
      final PriorBedRangeLoader p = new PriorBedRangeLoader();
      p.loadRanges(f);
      final ReferenceRanges<Double> referenceRanges = p.getReferenceRanges();
      final RangeList<Double> ranges = referenceRanges.get("chr1");
      final List<Double> v = ranges.find(1);
      assertEquals(1, v.size());
      assertEquals(0.42, v.get(0), 1e-12);
      assertEquals(1, ranges.find(99).size());
      assertNull(ranges.find(100));
    } finally {
      assertTrue(f.delete());
    }
  }
}
