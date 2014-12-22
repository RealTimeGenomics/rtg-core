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

import java.io.File;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.VcfReader;

import junit.framework.TestCase;

/**
 */
public class RocContainerTest extends TestCase {
  public void test() throws IOException {
    final File dir = FileUtils.createTempDir("rocCont", "test");
    try {
      final RocContainer roc = new RocContainer(RocSortOrder.DESCENDING);
      final File allFile = new File(dir, "all.txt.gz");
      roc.addFilter(RocFilter.ALL, allFile);
      final DetectedVariant v = new DetectedVariant(VcfReader.vcfLineToRecord(DetectedVariantTest.SHORT_LINE), 0, RocSortValueExtractor.NULL_EXTRACTOR, false);
      roc.addRocLine(0.1, 1.0, v);
      roc.addRocLine(0.1, 0.0, v);
      roc.addRocLine(0.2, 0.0, v);
      roc.addRocLine(0.1, 0.5, v);
      roc.addRocLine(0.2, 1.5, v);
      roc.addRocLine(0.1, 0.5, v);
      roc.addRocLine(0.3, 1.5, v);
      roc.writeRocs(5, true);
      final String all = FileHelper.gzFileToString(allFile);
      TestUtils.containsAll(all,
              "0.100\t5.000\t2",
              "0.200\t3.000\t1",
              "0.300\t1.500\t0");
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
}
