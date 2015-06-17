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

package com.rtg.simulation.reads;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class TaxonomyDistributionTest extends TestCase {

  @Override
  protected void tearDown() throws Exception {
    Diagnostic.setLogStream();
  }

  static final String EASY = ""
      + "# HEADER Line should be ignored" + StringUtils.LS
      + "0.5\t0" + StringUtils.LS
      + "0.5\t2" + StringUtils.LS
      ;
  public void testEasy() throws IOException {
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AAAAA", "AA", "AAA"));
    InputStream in = new ByteArrayInputStream(EASY.getBytes());
    Map<String, Integer> taxonLookup = new HashMap<>();
    taxonLookup.put("0", 0);
    taxonLookup.put("1", 1);
    taxonLookup.put("2", 2);
    MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    TaxonomyDistribution dist = new TaxonomyDistribution(in, taxonLookup, reader, TaxonomyDistribution.DistributionType.DNA_FRACTION);
    assertTrue(Arrays.equals(new double[]{0.5, 0, 0.5}, dist.getDistribution()));
    assertFalse(mps.toString().contains("Input distribution sums to"));

  }
  public void testHarder() throws IOException {
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AAAAA", "AA", "AAA", "AA"));
    InputStream in = new ByteArrayInputStream(EASY.getBytes());
    Map<String, Integer> taxonLookup = new HashMap<>();
    taxonLookup.put("0", 0);
    taxonLookup.put("1", 1);
    taxonLookup.put("2", 2);
    taxonLookup.put("3", 2);
    TaxonomyDistribution dist = new TaxonomyDistribution(in, taxonLookup, reader, TaxonomyDistribution.DistributionType.DNA_FRACTION);
    assertTrue(Arrays.equals(new double[]{0.5, 0, 0.3, 0.2}, dist.getDistribution()));

  }

  public void testParser() throws IOException {
    Map<Integer, Double> expected = new HashMap<>();
    expected.put(0, 0.5);
    expected.put(2, 0.5);
    assertEquals(expected, TaxonomyDistribution.parseTaxonDistribution(new ByteArrayInputStream(EASY.getBytes())));
  }
  public void testParserMalformed() {
    final String malformedLine = "0.5\t0\tI_don'tbelong";
    final String malformed = ""
                               + malformedLine  + StringUtils.LS
                               + "0.5\t2" + StringUtils.LS
        ;
    try {
      TaxonomyDistribution.parseTaxonDistribution(new ByteArrayInputStream(malformed.getBytes()));
      fail();
    } catch (IOException e) {
      assertEquals("Malformed line: " + malformedLine, e.getMessage());
    }
  }
  public void testParserNotDouble() {
    final String malformedLine = "0.monkey5\t0";
    final String malformed = ""
                             + malformedLine  + StringUtils.LS
                             + "0.5\t2" + StringUtils.LS
        ;
    try {
      TaxonomyDistribution.parseTaxonDistribution(new ByteArrayInputStream(malformed.getBytes()));
      fail();
    } catch (IOException e) {
      assertEquals("Malformed line: " + malformedLine, e.getMessage());
    }
  }
  public void testParserNotIntegral() {
    final String malformedLine = "0.5\t1.3";
    final String malformed = ""
                             + malformedLine  + StringUtils.LS
                             + "0.5\t2" + StringUtils.LS
        ;
    try {
      TaxonomyDistribution.parseTaxonDistribution(new ByteArrayInputStream(malformed.getBytes()));
      fail();
    } catch (IOException e) {
      assertEquals("Malformed line: " + malformedLine, e.getMessage());
    }
  }
  public void testParserDuplicateId() {
    final String malformedLine = "0.2\t2";
    final String malformed = ""
                             + "0.5\t2"  + StringUtils.LS
                             + malformedLine + StringUtils.LS
        ;
    try {
      TaxonomyDistribution.parseTaxonDistribution(new ByteArrayInputStream(malformed.getBytes()));
      fail();
    } catch (IOException e) {
      assertEquals("Duplicated key: " + malformedLine, e.getMessage());
    }
  }

  public void testDoesntSum() throws IOException {
    final String malformed = ""
                             + "0.49\t0"  + StringUtils.LS
                             + "0.49\t2" + StringUtils.LS
        ;
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AAAAA", "AA", "AAA"));
    InputStream in = new ByteArrayInputStream(malformed.getBytes());
    Map<String, Integer> taxonLookup = new HashMap<>();
    taxonLookup.put("0", 0);
    taxonLookup.put("1", 1);
    taxonLookup.put("2", 2);
    MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    TaxonomyDistribution dist = new TaxonomyDistribution(in, taxonLookup, reader, TaxonomyDistribution.DistributionType.DNA_FRACTION);
    assertTrue(Arrays.equals(new double[]{0.5, 0, 0.5}, dist.getDistribution()));
    TestUtils.containsAll(mps.toString(), "Input distribution sums to: 0.98");
  }

  public void testAbundance() throws IOException {
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AAAAAA", "AAA", "AA", "AAA"));
    InputStream in = new ByteArrayInputStream(EASY.getBytes());
    Map<String, Integer> taxonLookup = new HashMap<>();
    taxonLookup.put("0", 0);
    taxonLookup.put("1", 0);
    taxonLookup.put("2", 1);
    taxonLookup.put("3", 2);
    MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    TaxonomyDistribution dist = new TaxonomyDistribution(in, taxonLookup, reader, TaxonomyDistribution.DistributionType.ABUNDANCE);
    assertTrue(Arrays.toString(dist.getDistribution()), Arrays.equals(new double[]{0.5, 0.25, 0, 0.25}, dist.getDistribution()));
    assertFalse(mps.toString().contains("Input distribution sums to"));
  }
}
