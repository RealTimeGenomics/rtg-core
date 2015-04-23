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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 *
 */
public class VariantHelperTest extends TestCase {

  private static final String SNP = ""
    + VcfHeader.MINIMAL_HEADER + "\tFOO" + StringUtils.LS
    + "g1" + "\t" + "1" + "\t" + "." + "\t" + "A" + "\t" + "C" + "\t" + "3.0" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "2" + "\t" + "." + "\t" + "A" + "\t" + "G" + "\t" + "3.0" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "3" + "\t" + "." + "\t" + "A" + "\t" + "C" + "\t" + "3.8" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "4" + "\t" + "." + "\t" + "T" + "\t" + "G" + "\t" + "3.8" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "4" + "\t" + "." + "\t" + "A" + "\t" + "C" + "\t" + "4.7" + "\t" + "PASS\t.\tGT\t0/0" + StringUtils.LS
    + "g1" + "\t" + "5" + "\t" + "." + "\t" + "C" + "\t" + "G" + "\t" + "3.9" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "6" + "\t" + "." + "\t" + "G" + "\t" + "C" + "\t" + "4.5" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
    + "g1" + "\t" + "7" + "\t" + "." + "\t" + "A" + "\t" + "G" + "\t" + "4.7" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
   ;

  public final void testLoadSnpRange() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("loadsnps", "commonserverutils");
    try {
      final File snps = new File(f, "snps");
      FileUtils.stringToFile(SNP, snps);
      checkResults(snps);
      final File snpsgz = new File(f, "snps.gz");
      FileHelper.stringToGzFile(SNP, snpsgz);
      checkResults(snpsgz);
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  private void checkResults(final File snps) throws IOException {
    final LinkedHashMap<String, ArrayList<AviewVariant>> snpmap = new LinkedHashMap<>();
    VariantHelper.loadSnpRange(snpmap, snps,  new RegionRestriction("g1", 1, 4));
    //this should be 3 as the fourth variant at position 4 is same as reference and should be discarded
    assertEquals(snpmap.get("FOO").toString(), 3, snpmap.get("FOO").size());
  }

  private static final String SNP_NO_SAMPLE = ""
      + VcfHeader.MINIMAL_HEADER + StringUtils.LS
      + "g1" + "\t" + "1" + "\t" + "." + "\t" + "A" + "\t" + "T" + "\t" + "3.0" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
      + "g1" + "\t" + "7" + "\t" + "." + "\t" + "A" + "\t" + "T" + "\t" + "4.7" + "\t" + "PASS\t.\tGT\t0/1" + StringUtils.LS
     ;

  public void testNoSampleVCF() throws IOException {
    final File f = FileUtils.createTempDir("loadsnps", "commonserverutils");
    try {
      final File snps = new File(f, "snps");
      FileUtils.stringToFile(SNP_NO_SAMPLE, snps);
      try {
        final LinkedHashMap<String, ArrayList<AviewVariant>> snpmap = new LinkedHashMap<>();
        VariantHelper.loadSnpRange(snpmap, snps, new RegionRestriction("g1", 1, 4));
        fail();
      } catch (NoTalkbackSlimException ex) {
        assertEquals("VCF header line contains format field but no sample fields", ex.getMessage());
      }

    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }

  }

  private static final String SNP_MULTI_SAMPLE = ""
      + VcfHeader.MINIMAL_HEADER + "\t" + "SAMPLE1" + "\t" + "SAMPLE2" + StringUtils.LS
      + "g1" + "\t" + "3" + "\t" + "." + "\t" + "A" + "\t" + "T" + "\t" + "3.0" + "\t" + "PASS\t.\tGT\t0/1\t0/1" + StringUtils.LS
      + "g1" + "\t" + "7" + "\t" + "." + "\t" + "A" + "\t" + "T" + "\t" + "4.7" + "\t" + "PASS\t.\tGT\t0/1\t1/1" + StringUtils.LS
     ;

  public void testMultiSampleVCF() throws IOException {
    final File f = FileUtils.createTempDir("loadsnps", "commonserverutils");
    try {
      final File snps = new File(f, "snps");
      FileUtils.stringToFile(SNP_MULTI_SAMPLE, snps);
      try {
        final LinkedHashMap<String, ArrayList<AviewVariant>> snpmap = new LinkedHashMap<>();
        VariantHelper.loadSnpRange(snpmap, snps,  new RegionRestriction("g1", 1, 4));
        fail();
      } catch (NoTalkbackSlimException ex) {
        assertEquals("No sample name specified but multiple samples available. Please select from the samples available: SAMPLE1 SAMPLE2", ex.getMessage());
      }

    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }

  }

  public void testRetriveMultiSampleVCF() throws IOException {
    final File f = FileUtils.createTempDir("loadsnps", "commonserverutils");
    try {
      final File snps = new File(f, "snps");
      FileUtils.stringToFile(SNP_MULTI_SAMPLE, snps);
      final LinkedHashMap<String, ArrayList<AviewVariant>> snpmap = new LinkedHashMap<>();
      VariantHelper.loadSnpRange(snpmap, snps,  new RegionRestriction("g1", 1, 4), "SAMPLE1", "SAMPLE2");
      assertEquals(snpmap.get("SAMPLE1").toString(), 1, snpmap.get("SAMPLE1").size());
      assertEquals(snpmap.get("SAMPLE2").toString(), 1, snpmap.get("SAMPLE2").size());
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }

  }
}
