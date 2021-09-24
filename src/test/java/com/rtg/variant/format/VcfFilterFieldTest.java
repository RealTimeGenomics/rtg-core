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

package com.rtg.variant.format;

import java.io.File;

import com.rtg.util.TestUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfFilterFieldTest extends TestCase {

  public void testEnum() {
    TestUtils.testEnum(VcfFilterField.class, "[OC, A, RC, RX, RCEQUIV, IONT, BED, PASS, OTHER]");
    assertEquals(VcfFilterField.values().length - 1, VcfFilterField.OTHER.ordinal());
  }

  public void testHeaders() {
    final VariantParams params = VariantParams.builder().maxAmbiguity(0.1).regionsFilterBedFile(new File("foo.bed")).create();
    final VcfHeader header = new VcfHeader();
    for (VcfFilterField field : VcfFilterField.values()) {
      field.updateHeader(header, params);
    }
    final String expected = ""
    + "##FILTER=<ID=OC,Description=\"Coverage threshold exceeded\">\n"
    + "##FILTER=<ID=a10.0,Description=\"Ambiguity exceeded 10.0\">\n"
    + "##FILTER=<ID=RC,Description=\"RTG variant is a complex region\">\n"
    + "##FILTER=<ID=RX,Description=\"RTG variant contained within hypercomplex region\">\n"
    + "##FILTER=<ID=RCEQUIV,Description=\"RTG variant is equivalent to the previous variant\">\n"
    + "##FILTER=<ID=IONT,Description=\"IonTorrent specific filter applied\">\n"
    + "##FILTER=<ID=BED,Description=\"Variant falls outside the specified target regions\">\n"
    + "##FILTER=<ID=OTHER,Description=\"Variant is invalid for unknown reasons\">\n"
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    assertEquals(expected, header.toString());
  }

  public void testDummyAmbiguity() {
    final VcfHeader header = new VcfHeader();
    VcfFilterField.A.updateHeader(header, null);
    final String expected = ""
    + "##FILTER=<ID=a-1.0,Description=\"Ambiguity exceeded -1.0\">\n"
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    assertEquals(expected, header.toString());
  }

  public void testFilters() {
    final VariantParams params = VariantParams.builder().maxAmbiguity(0.1).create();
    final Variant call = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    for (VariantFilter filter : VariantFilter.values()) {
      call.addFilter(filter);
    }
    final VcfRecord rec = new VcfRecord("foo", 0, "c");
    for (VcfFilterField field : VcfFilterField.values()) {
      field.updateRecord(rec, call, params);
    }

    assertEquals("foo\t1\t.\tc\t.\t.\tOC;a10.0;RX;RCEQUIV;IONT;BED\t.", rec.toString());
  }

  public void testComplexRegionFilter() {
    final VariantParams params = VariantParams.builder().create();
    final Variant call = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    call.addFilter(VariantFilter.FAILED_COMPLEX);
    final VcfRecord rec = new VcfRecord("foo", 0, "c");
    VcfFilterField.RC.updateRecord(rec, call, params);

    assertEquals("foo\t1\t.\tc\t.\t.\tRC\t.", rec.toString());
  }

  public void testOtherFilter() {
    final VariantParams params = VariantParams.builder().create();
    final Variant call = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    call.addFilter(VariantFilter.OTHER);
    final VcfRecord rec = new VcfRecord("foo", 0, "c");
    VcfFilterField.OTHER.updateRecord(rec, call, params);

    assertEquals("foo\t1\t.\tc\t.\t.\tOTHER\t.", rec.toString());
  }

  public void testPassFilter() {
    final VariantParams params = VariantParams.builder().create();
    final Variant call = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    final VcfRecord rec = new VcfRecord("foo", 0, "c");
    VcfFilterField.PASS.updateRecord(rec, call, params);

    assertEquals("foo\t1\t.\tc\t.\t.\tPASS\t.", rec.toString());
  }
}
