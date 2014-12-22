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
import static com.rtg.util.StringUtils.TAB;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Collections;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfReaderTest extends TestCase {

  public void testEmpty() {
    checkIoException("", "file format");
  }

  static final String[] BAD_HDR = {
    "#CHRM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "#CHROM\tPOS\tID\tREF\tALT\tqual\tFILTER\tINFO",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", // needs some sample IDs
  };

  public void testBadHdr() throws Exception {
    for (final String h : BAD_HDR) {
      final Reader in = new StringReader(""
        + VcfHeader.VERSION_LINE + LS
        + h + LS);
      try {
        new VcfReader(new BufferedReader(in)).hasNext();
        fail();
      } catch (final NoTalkbackSlimException ex) {
        TestUtils.containsAll(ex.getMessage(), "VCF", "Header");
      }
    }
  }

  protected void checkIoException(String header, String... contains) {
    final Reader in = new StringReader(header);
    try {
      new VcfReader(new BufferedReader(in)).hasNext();
      fail("IOException expected from .vcf file: " + header);
    } catch (final IOException ex) {
      TestUtils.containsAll(ex.getMessage(), contains);
    }
  }

  static final String[] BAD_RECORD = {
    "chr1   123    .    G    A    29    PASS",                                 "record",
    "chr1   123    foo  G    A    29    PASS    X=yy;DP=7    GT:GQ",           "record",
    "chr1   123    foo  G    A    29    PASS    NS=3;DP=7    GT:GQ   0|0:34",  "record",
  };

  public void testBadRecords() {
    for (int i = 0; i < BAD_RECORD.length; i += 2) {
      final String badrec = BAD_RECORD[i].replaceAll("  *", TAB);
      checkIoException(HEADER0 + badrec + LS,
        "Invalid VCF record");
    }
  }

  /** standard header to share */
  public static final String HEADER0 = ""
    + "##fileformat=VCFv4.1" + LS
    + "#CHROM" + TAB + "POS" + TAB + "ID" + TAB + "REF" + TAB + "ALT" + TAB + "QUAL" + TAB + "FILTER" + TAB
        + "INFO" + LS
    ;
  /** standard header to share */
  public static final String HEADER0_B = ""
    + "##fileformat=VCFv4.1" + LS
    + "#CHROM" + TAB + "POS" + TAB + "ID" + TAB + "REF" + TAB + "ALT" + TAB + "QUAL" + TAB + "FILTER" + TAB
        + "INFO" + TAB + "FORMAT" + TAB + "SAMPLE" + LS
    ;
  public void testHeader0() throws IOException {
    final Reader in = new StringReader(HEADER0);
    final VcfReader reader = new VcfReader(new BufferedReader(in));
    final VcfHeader hdr = reader.getHeader();
    assertNotNull(hdr);
    assertEquals(0, hdr.getGenericMetaInformationLines().size());
    assertEquals(0, hdr.getNumberOfSamples());
    assertFalse(reader.hasNext());
  }

  protected static final String HEADER1 = ""
    + "##fileformat=VCFv4.1" + LS
    + "##X=YYY" + LS
    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" + LS
    + "#CHROM" + TAB + "POS" + TAB + "ID" + TAB + "REF" + TAB + "ALT" + TAB + "QUAL" + TAB + "FILTER" + TAB + "INFO" + TAB + "FORMAT" + TAB + "sample1" + LS
    + "chr1" + TAB + "123" + TAB + "." + TAB + "G" + TAB + "A,G" + TAB + "29" + TAB + "q3;s5;c" + TAB + "X=yy;DP=7" + TAB + "GT:GQ" + TAB + "0|0:7" + LS
    ;
  public void testNext() throws IOException {
    final Reader in = new StringReader(HEADER1);
    final VcfReader reader = new VcfReader(new BufferedReader(in));
    final VcfHeader hdr = reader.getHeader();
    assertNotNull(hdr);
    assertEquals(2, hdr.getGenericMetaInformationLines().size() + hdr.getInfoLines().size() + hdr.getFormatLines().size() + hdr.getFilterLines().size());
    assertEquals(1, hdr.getNumberOfSamples());
    assertTrue(reader.hasNext());
    final VcfRecord rec = reader.next();
    assertNotNull(rec);
    assertEquals("chr1", rec.getSequenceName());
    assertEquals(123, rec.getOneBasedStart());
    assertEquals(".", rec.getId());
    assertEquals("G", rec.getRefCall());

    assertEquals(2, rec.getAltCalls().size());
    assertEquals("A", rec.getAltCalls().get(0));
    assertEquals("G", rec.getAltCalls().get(1));
    assertEquals("29", rec.getQuality());

    assertEquals(3, rec.getFilters().size());
    assertEquals("q3", rec.getFilters().get(0));
    assertEquals("s5", rec.getFilters().get(1));
    assertEquals("c", rec.getFilters().get(2));

    assertEquals(2, rec.getInfo().size());
    assertEquals("yy", rec.getInfo().get("X").iterator().next());
    assertEquals("7", rec.getInfo().get("DP").iterator().next());

    assertEquals(2, rec.getFormatAndSample().size());
    assertEquals(Collections.singletonList("0|0"), rec.getFormatAndSample().get("GT"));
    assertEquals(Collections.singletonList("7"), rec.getFormatAndSample().get("GQ"));
  }
}
