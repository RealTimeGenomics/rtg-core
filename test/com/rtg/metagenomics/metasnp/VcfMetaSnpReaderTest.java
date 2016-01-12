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
package com.rtg.metagenomics.metasnp;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class VcfMetaSnpReaderTest extends TestCase {

  private static final String HEADER = ""
    + "##fileformat=VCFv4.1" + StringUtils.LS
    + "##fileDate=20151211" + StringUtils.LS
    + "##contig=<ID=chr1,length=86738>" + StringUtils.LS
    + "##FORMAT=<ID=ADE,Number=.,Type=Float,Description=\"Allelic depths\">" + StringUtils.LS;

  private static final String SINGLE_SAMPLE_VCF_FILE = HEADER
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor" + StringUtils.LS
    + "chr1\t3\t.\tC\tG\t84.4\tPASS\t.\tAD\t1.9,6.0" + StringUtils.LS
    + "chr1\t239\t.\tA\tT\t0.0\tPASS\t.\tAD\t250.3,231.2" + StringUtils.LS;

  public void testSingleSample() throws IOException {
    try (VcfMetaSnpReader reader = new VcfMetaSnpReader(new BufferedReader(new StringReader(SINGLE_SAMPLE_VCF_FILE)))) {
      final MetaSnpLine line = reader.nextLine();
      assertEquals("chr1", line.getSequence());
      assertEquals(2, line.getPosition());
      assertEquals(0, line.getReferenceIndex());
      assertEquals("C", line.getReferenceAllele());
      final double[][] cnts = line.mCounts;
      assertEquals(2, cnts.length);
      assertEquals(1.9, cnts[0][0], 1e-9);
      assertEquals(6.0, cnts[1][0], 1e-9);
      final MetaSnpLine line2 = reader.nextLine();
      assertEquals("chr1", line2.getSequence());
      assertEquals(238, line2.getPosition());
      assertEquals(0, line2.getReferenceIndex());
      assertEquals("A", line2.getReferenceAllele());
      final double[][] cnts1 = line2.mCounts;
      assertEquals(2, cnts1.length);
      assertEquals(250.3, cnts1[0][0], 1e-9);
      assertEquals(231.2, cnts1[1][0], 1e-9);
      assertNull(reader.nextLine());
    }
  }

  private static final String TWO_SAMPLE_VCF_FILE = HEADER
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\tcancer" + StringUtils.LS
    + "chr1\t3\t.\tC\tG\t84.4\tPASS\t.\tAD\t1.9,6.0\t3.1,4.1" + StringUtils.LS
    + "chr1\t239\t.\tA\tT\t0.0\tPASS\t.\tAD\t250.3,231.2\t7,42" + StringUtils.LS;

  public void testTwoSample() throws IOException {
    try (VcfMetaSnpReader reader = new VcfMetaSnpReader(new BufferedReader(new StringReader(TWO_SAMPLE_VCF_FILE)))) {
      final MetaSnpLine line = reader.nextLine();
      assertEquals("chr1", line.getSequence());
      assertEquals(2, line.getPosition());
      assertEquals(0, line.getReferenceIndex());
      assertEquals("C", line.getReferenceAllele());
      final double[][] cnts = line.mCounts;
      assertEquals(2, cnts.length);
      assertEquals(1.9, cnts[0][0], 1e-9);
      assertEquals(6.0, cnts[1][0], 1e-9);
      assertEquals(3.1, cnts[0][1], 1e-9);
      assertEquals(4.1, cnts[1][1], 1e-9);
      final MetaSnpLine line2 = reader.nextLine();
      assertEquals("chr1", line2.getSequence());
      assertEquals(238, line2.getPosition());
      assertEquals(0, line2.getReferenceIndex());
      assertEquals("A", line2.getReferenceAllele());
      final double[][] cnts1 = line2.mCounts;
      assertEquals(2, cnts1.length);
      assertEquals(250.3, cnts1[0][0], 1e-9);
      assertEquals(231.2, cnts1[1][0], 1e-9);
      assertEquals(7, cnts1[0][1], 1e-9);
      assertEquals(42, cnts1[1][1], 1e-9);
      assertNull(reader.nextLine());
    }
  }

  private static final String MISSING_AD_VCF_FILE = HEADER
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor" + StringUtils.LS
    + "chr1\t3\t.\tC\tG\t84.4\tPASS\t.\tAD\t." + StringUtils.LS
    + "chr1\t239\t.\tA\tT\t0.0\tPASS\t.\tAD\t250.3,231.2" + StringUtils.LS;

  public void testMissingAdSample() throws IOException {
    try (VcfMetaSnpReader reader = new VcfMetaSnpReader(new BufferedReader(new StringReader(MISSING_AD_VCF_FILE)))) {
      final MetaSnpLine line2 = reader.nextLine();
      final double[][] cnts1 = line2.mCounts;
      assertEquals(2, cnts1.length);
      assertEquals(250.3, cnts1[0][0], 1e-9);
      assertEquals(231.2, cnts1[1][0], 1e-9);
      assertNull(reader.nextLine());
    }
  }

}
