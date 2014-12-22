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
import java.util.Set;

import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfFormatStripperTest extends TestCase {

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
            .addInfo("AN", "5")
            .setNumberOfSamples(2)
            .addFormatAndSample("GT", "0/0")
            .addFormatAndSample("GT", "0/1")
            .addFormatAndSample("DS", "je")
            .addFormatAndSample("DS", "er")
            .addFormatAndSample("GL", "je")
            .addFormatAndSample("GL", "fe")
            .addFormatAndSample("GQ", "100")
            .addFormatAndSample("GQ", "95")
    ;
    return rec;
  }

  private VcfHeader createTestHeader() {
    final VcfHeader head = new VcfHeader();
    head.addFormatField(VcfHeader.parseFormatLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    head.addFormatField(VcfHeader.parseFormatLine("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage from MaCH/Thunder\">"));
    head.addFormatField(VcfHeader.parseFormatLine("##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihoods\">"));
    return head;
  }

  public void testKeep() {
    final HashSet<String> formats = new HashSet<>();
    formats.add("GT");

    final VcfFormatStripper ann = new VcfFormatStripper(formats, true);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final Set<String> formatlist = rec.getFormats();
    assertNotNull(formatlist);
    assertEquals(1, formatlist.size());
    assertTrue(formatlist.contains("GT"));

    final Map<String, ArrayList<String>> formatMap = rec.getFormatAndSample();
    assertNotNull(formatMap);
    assertEquals(1, formatMap.size());

    assertTrue(formatMap.containsKey("GT"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<FormatField> headerformats = header.getFormatLines();
    assertNotNull(headerformats);
    assertEquals(1, headerformats.size());

    assertEquals("GT", headerformats.get(0).getId());
  }

  public void testRemove() {
    final HashSet<String> formats = new HashSet<>();
    formats.add("DS");
    formats.add("GL");

    final VcfFormatStripper ann = new VcfFormatStripper(formats, false);

    final VcfRecord rec = createTestRecord();
    ann.annotate(rec);

    final Set<String> formatset = rec.getFormats();
    assertNotNull(formatset);

    assertEquals(2, formatset.size());
    assertTrue(formatset.contains("GT"));
    assertTrue(formatset.contains("GQ"));

    final Map<String, ArrayList<String>> formatMap = rec.getFormatAndSample();
    assertNotNull(formatMap);
    assertEquals(2, formatMap.size());

    assertTrue(formatMap.containsKey("GT"));
    assertTrue(formatMap.containsKey("GQ"));

    final VcfHeader header = createTestHeader();
    ann.updateHeader(header);

    final List<FormatField> formatLines = header.getFormatLines();
    assertNotNull(formatLines);
    assertEquals(1, formatLines.size()); //didn't add GQ in, so only 1

    assertEquals("GT", formatLines.get(0).getId());
  }
}
