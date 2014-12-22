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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfInfoStripperTest extends TestCase {

  private VcfRecord createTestRecord() {
    final VcfRecord rec = new VcfRecord();
    rec.setSequence("chr1")
            .setStart(1209)
            .setId(".")
            .setQuality("12.8")
            .setRefCall("a")
            .addAltCall("c")
            .addAltCall("t")
            .addFilter("TEST1")
            .addFilter("TEST2")
            .addInfo("DP", "23")
            .addInfo("AC", "2")
            .addInfo("TEST", "45", "46", "47", "48")
            .addInfo("AN", "5")
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
    head.addInfoField(VcfHeader.parseInfoLine("##INFO=<ID=yo, Number=5, Type=Float, Description=\"fun for the whole family\">"));
    head.addInfoField(VcfHeader.parseInfoLine("##INFO=<ID=no, Number=5, Type=Float, Description=\"fun for the whole family\">"));
    head.addInfoField(VcfHeader.parseInfoLine("##INFO=<ID=AC, Number=5, Type=Float, Description=\"fun for the whole family\">"));
    head.addInfoField(VcfHeader.parseInfoLine("##INFO=<ID=go, Number=5, Type=Float, Description=\"fun for the whole family\">"));
    head.addInfoField(VcfHeader.parseInfoLine("##INFO=<ID=AN, Number=5, Type=Float, Description=\"fun for the whole family\">"));
    return head;
  }

  public void testKeep() {
    final HashSet<String> infos = new HashSet<>();
    infos.add("AN");
    infos.add("AC");

    final VcfInfoStripper ann = new VcfInfoStripper(infos, true);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final Map<String, ArrayList<String>> infoMap = rec.getInfo();
    assertNotNull(infoMap);
    assertEquals(2, infoMap.size());

    assertTrue(infoMap.containsKey("AN"));
    assertTrue(infoMap.containsKey("AC"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<InfoField> headerinfos = header.getInfoLines();
    assertNotNull(headerinfos);
    assertEquals(2, headerinfos.size());

    assertEquals("AC", headerinfos.get(0).getId());
    assertEquals("AN", headerinfos.get(1).getId());
  }

  public void testRemove() {
    final HashSet<String> infos = new HashSet<>();
    infos.add("AN");
    infos.add("AC");

    final VcfInfoStripper ann = new VcfInfoStripper(infos, false);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final Map<String, ArrayList<String>> infoMap = rec.getInfo();
    assertNotNull(infoMap);
    assertEquals(2, infoMap.size());

    assertTrue(infoMap.containsKey("DP"));
    assertTrue(infoMap.containsKey("TEST"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<InfoField> headerinfos = header.getInfoLines();
    assertNotNull(headerinfos);
    assertEquals(3, headerinfos.size());

    assertEquals("yo", headerinfos.get(0).getId());
    assertEquals("no", headerinfos.get(1).getId());
    assertEquals("go", headerinfos.get(2).getId());
  }

  public void testRemoveAll() {
    final VcfInfoStripper ann = new VcfInfoStripper(true);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final Map<String, ArrayList<String>> infoMap = rec.getInfo();
    assertNotNull(infoMap);
    assertEquals(0, infoMap.size());

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<InfoField> infos = header.getInfoLines();
    assertNotNull(infos);
    assertEquals(0, infos.size());
  }
}
