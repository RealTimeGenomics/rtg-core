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

import java.util.HashSet;
import java.util.List;

import com.rtg.vcf.header.FilterField;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfFilterStripperTest extends TestCase {

  private VcfRecord createTestRecord() {
    final VcfRecord rec = new VcfRecord();
    rec.setSequence("chr1")
            .setStart(1209)
            .setId(".")
            .setQuality("12.8")
            .setRefCall("a")
            .addAltCall("c")
            .addAltCall("t")
            .addFilter("no")
            .addFilter("wo")
            .addFilter("go")
            .addInfo("DP", "23")
            .addInfo("TEST", "45", "46", "47", "48")
            .setNumberOfSamples(2)
            .addFormatAndSample("GT", "0/0")
            .addFormatAndSample("GT", "0/1")
            .addFormatAndSample("GQ", "100")
            .addFormatAndSample("GQ", "95")
    ;
    return rec;
  }

  private VcfHeader createTestHeader() {
    final VcfHeader head = new VcfHeader();
    head.addFilterField("no", "blah");
    head.addFilterField("wo", "blah");
    head.addFilterField("go", "blah");
    return head;
  }

  public void testKeep() {
    final HashSet<String> filters = new HashSet<>();
    filters.add("wo");

    final VcfFilterStripper ann = new VcfFilterStripper(filters, true);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final List<String> filterList = rec.getFilters();
    assertNotNull(filterList);
    assertEquals(1, filterList.size());

    assertTrue(filterList.get(0).equals("wo"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<FilterField> headerfilters = header.getFilterLines();
    assertNotNull(headerfilters);
    assertEquals(1, headerfilters.size());

    assertEquals("wo", headerfilters.get(0).getId());
  }

  public void testRemove() {
    final HashSet<String> filters = new HashSet<>();
    filters.add("no");
    filters.add("go");

    final VcfFilterStripper ann = new VcfFilterStripper(filters, false);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final List<String> filterlist = rec.getFilters();
    assertNotNull(filterlist);
    assertEquals(1, filterlist.size());

    assertTrue(filterlist.get(0).equals("wo"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<FilterField> headerfilters = header.getFilterLines();
    assertNotNull(headerfilters);
    assertEquals(1, headerfilters.size());

    assertEquals("wo", headerfilters.get(0).getId());
  }

  public void testRemoveAll() {
    final VcfFilterStripper ann = new VcfFilterStripper(true);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final List<String> filterlist = rec.getFilters();
    assertNotNull(filterlist);
    assertEquals(0, filterlist.size());

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<FilterField> headerfilters = header.getFilterLines();
    assertNotNull(headerfilters);
    assertEquals(0, headerfilters.size());
  }
}
