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

package com.rtg.variant.bayes.multisample;

import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VariantOutputVcfFormatterTest;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class VcfEquivalentFilterTest extends TestCase {

  private Variant getVariant(final int start, final int end, final String bestName) {
    final StringBuilder sb = new StringBuilder();
    for (int i = start; i < end; ++i) {
      sb.append("A");
    }
    final String refNts = sb.toString();
    final VariantSample vs = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, bestName, refNts.equals(bestName), 40.0, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final VariantLocus locus = new VariantLocus("foo", start, end, refNts, '-');
    final Variant v = new Variant(locus, vs);
    v.setInteresting();
    v.setComplexScored();
    return v;
  }

  private VcfRecord getCalls(int start, int end, String bestName) {
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter("foo");
    return formatter.makeVcfRecord(getVariant(start, end, bestName));
  }

  private VcfRecord getMultiSampleCalls(int start, int end, String bestName, String secondName) {
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter("foo", "bar");
    final Variant subvar1 = getVariant(start, end, bestName);
    final VariantLocus locus = new VariantLocus("foo", start, end, subvar1.getLocus().getRefNts(), '-');
    final Variant v = new Variant(locus, subvar1.getSample(0), secondName == null ? null : getVariant(start, end, secondName).getSample(0));
    v.setInteresting();
    v.setComplexScored();
    final VcfRecord rec = formatter.makeVcfRecord(v);
    assertFalse(isComplexEquivalent(rec));
    return rec;
  }

  private boolean isComplexEquivalent(final VcfRecord rec) {
    return rec.hasInfo("RCE");
  }

  public void testSimpleHomozygousAlt() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(2, 3, "AA"));
    list.add(getCalls(7, 8, "AA"));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(isComplexEquivalent(list.get(1)));
  }

  public void testSimpleHomozygousRef() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(2, 3, "A"));
    list.add(getCalls(7, 8, "A"));
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSimpleHeterozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(2, 3, "A:AA"));
    list.add(getCalls(7, 8, "A:AA"));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(isComplexEquivalent(list.get(1)));
  }

  public void testSmallerSeparationHomozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(2, 3, "AA"));
    list.add(getCalls(7, 8, "AA"));
    f.filter(list, 3);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSmallerSeparationHheterozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(2, 3, "A:AA"));
    list.add(getCalls(7, 8, "A:AA"));
    f.filter(list, 3);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSimpleHomozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "AA"));
    list.add(getMultiSampleCalls(7, 8, "AA", "AA"));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(isComplexEquivalent(list.get(1)));
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "AA", "AA"));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(list.size() == 1);
  }

  public void testSimpleHeterozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "A:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "A:AA"));

    final List<VcfRecord> filteredCalls = f.filter(list, 10); //without last one
    assertTrue(isComplexEquivalent(filteredCalls.get(0)));
    assertTrue(isComplexEquivalent(f.lastCall().get(0)));
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "A:AA", "A:AA"));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
  }

  public void testSmallerSeparationHomozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "AA"));
    list.add(getMultiSampleCalls(7, 8, "AA", "AA"));
    f.filter(list, 3);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSmallerSeparationHeterozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "A:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "A:AA"));
    f.filter(list, 3);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));

  }

  public void testSimpleHomozygousMultisampleDifferent() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "TT"));
    list.add(getMultiSampleCalls(7, 8, "AA", "GG"));
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSimpleHeterozygousMultisampleDifferent() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "T:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "T:TT"));
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testSimpleHomozygousFirstNull() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(7, 8, "AA"));
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
  }

  public void testSimpleHomozygousNameNull() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getCalls(7, 8, "A"));
    list.add(getCalls(10, 11, "T"));
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
  }

  public void testMultisampleYChr() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", null));
    list.add(getMultiSampleCalls(7, 8, "AA", null));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(isComplexEquivalent(list.get(1)));
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "AA", null));
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(list.size() == 1);
  }

  public void testMultisampleParCross() { // Unlikely unless using an odd sex chromosome configuration
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final VcfEquivalentFilter f = new VcfEquivalentFilter(template, null);
    final List<VcfRecord> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "AA")); // Male in non-PAR Y
    list.add(getMultiSampleCalls(7, 8, "AA", null)); // Male inside PAR Y
    f.filter(list, 10);
    assertFalse(isComplexEquivalent(list.get(0)));
    assertFalse(isComplexEquivalent(list.get(1)));
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "AA", null)); // Male inside PAR Y
    f.filter(list, 10);
    assertTrue(isComplexEquivalent(list.get(0)));
    assertTrue(list.size() == 1);
  }

}
