/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
