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

package com.rtg.tabix;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.sam.SamRangeUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 * Test class
 */
public class TabixLineReaderTest extends TestCase {

  protected NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    Diagnostic.setLogStream();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public static String extractRecords(File input, File tabix, RegionRestriction range) throws IOException {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final OutputStream out = mps.outputStream();
      try (TabixLineReader reader = new TabixLineReader(input, tabix, range)) {
        String line;
        while ((line = reader.readLine()) != null) {
          out.write(line.getBytes());
          out.write(StringUtils.LS.getBytes());
        }
      }
      return mps.toString();
    }
  }

  public static String extractRecords(File input, File tabix, ReferenceRanges<String> ranges) throws IOException {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final OutputStream out = mps.outputStream();
      try (TabixLineReader reader = new TabixLineReader(input, tabix, ranges)) {
        String line;
        while ((line = reader.readLine()) != null) {
          out.write(line.getBytes());
          out.write(StringUtils.LS.getBytes());
        }
      }
      return mps.toString();
    }
  }

  public void testSingleRegion() throws IOException {
    try (final TestDirectory dir = new TestDirectory("TabixLineReader")) {

      final File input = new File(dir, "snp_only.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input);
      final File tabix = new File(dir, "snp_only.vcf.gz.tbi");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz.tbi", tabix);

      final RegionRestriction region = new RegionRestriction("simulatedSequence19", 500, 1000);
      String result = extractRecords(input, tabix, region);
      mNano.check("tlr-single-region", result);

      final ReferenceRanges<String> ranges = SamRangeUtils.createExplicitReferenceRange(region);
      String result2 = extractRecords(input, tabix, ranges);
      assertEquals(result, result2);
    }
  }

  public void testMultiRegion() throws IOException {
    try (final TestDirectory dir = new TestDirectory("TabixLineReader")) {

      final File input = new File(dir, "snp_only.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input);
      final File tabix = new File(dir, "snp_only.vcf.gz.tbi");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz.tbi", tabix);
      ReferenceRanges<String> ranges = SamRangeUtils.createExplicitReferenceRange(
        new RegionRestriction("simulatedSequence2", 215, 3345),
        new RegionRestriction("simulatedSequence14", 0, 1567),
        new RegionRestriction("simulatedSequence14", 1567, 10000),
        new RegionRestriction("simulatedSequence19")
        );
      String result = extractRecords(input, tabix, ranges);
      mNano.check("tlr-multi-region", result);
    }
  }
}
