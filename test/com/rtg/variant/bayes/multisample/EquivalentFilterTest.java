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

package com.rtg.variant.bayes.multisample;

import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.format.VariantOutputVcfFormatterTest;

import junit.framework.TestCase;

/**
 */
public class EquivalentFilterTest extends TestCase {

  private Variant getCalls(int start, int end, String bestName) {
    StringBuilder sb = new StringBuilder();
    for (int i = start; i < end; i++) {
      sb.append("A");
    }
    final VariantSample vs = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, bestName, false, 40.0, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final VariantLocus locus = new VariantLocus("foo", start, end, sb.toString(), '-');
    final Variant v = new Variant(locus, vs);
    v.setInteresting();
    v.setComplexScored();
    return v;
  }

  private Variant getMultiSampleCalls(int start, int end, String bestName, String secondName) {
    final Variant subvar1 = getCalls(start, end, bestName);
    final Variant subvar2 = getCalls(start, end, secondName);
    final VariantLocus locus = new VariantLocus("foo", start, end, subvar1.getLocus().getRefNts(), '-');
    final Variant v = new Variant(locus, subvar1.getSample(0), subvar2.getSample(0));
    v.setInteresting();
    v.setComplexScored();
    return v;
  }

  public void testSimpleHomozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(2, 3, "AA"));
    list.add(getCalls(7, 8, "AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertTrue(list.get(0).isComplexEquivalent());
    assertTrue(list.get(1).isComplexEquivalent());
  }

  public void testSimpleHeterozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(2, 3, "A:AA"));
    list.add(getCalls(7, 8, "A:AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertTrue(list.get(0).isComplexEquivalent());
    assertTrue(list.get(1).isComplexEquivalent());
  }

  public void testSmallerSeparationHomozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(2, 3, "AA"));
    list.add(getCalls(7, 8, "AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 3);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }

  public void testSmallerSeparationHheterozygous() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(2, 3, "A:AA"));
    list.add(getCalls(7, 8, "A:AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 3);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }

  public void testSimpleHomozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "AA"));
    list.add(getMultiSampleCalls(7, 8, "AA", "AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertTrue(list.get(0).isComplexEquivalent());
    assertTrue(list.get(1).isComplexEquivalent());
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "AA", "AA"));
    f.filter(list, 10);
    assertTrue(list.get(0).isComplexEquivalent());
    assertTrue(list.size() == 1);
  }

  public void testSimpleHeterozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "A:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "A:AA"));

    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());

    final List<Variant> filteredCalls = f.filter(list, 10); //without last one
    assertTrue(filteredCalls.get(0).isComplexEquivalent());
    assertTrue(f.lastCall().get(0).isComplexEquivalent());
    list.clear();
    list.add(getMultiSampleCalls(12, 13, "A:AA", "A:AA"));
    f.filter(list, 10);
    assertTrue(list.get(0).isComplexEquivalent());
  }

  public void testSmallerSeparationHomozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "AA"));
    list.add(getMultiSampleCalls(7, 8, "AA", "AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 3);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }

  public void testSmallerSeparationHeterozygousMultisample() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "A:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "A:AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 3);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());

  }

  public void testSimpleHomozygousMultisampleDifferent() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "AA", "TT"));
    list.add(getMultiSampleCalls(7, 8, "AA", "GG"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }

  public void testSimpleHeterozygousMultisampleDifferent() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getMultiSampleCalls(2, 3, "A:AA", "T:AA"));
    list.add(getMultiSampleCalls(7, 8, "A:AA", "T:TT"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }

  public void testSimpleHomozygousFirstNull() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(7, 8, "AA"));
    assertFalse(list.get(0).isComplexEquivalent());
    f.filter(list, 10);
    assertFalse(list.get(0).isComplexEquivalent());
  }

  public void testSimpleHomozygousNameNull() {
    final byte[] template = DNA.stringDNAtoByte("AAAAAAAAAAAAAAAAAAAA");
    final EquivalentFilter f = new EquivalentFilter(template, null);
    final List<Variant> list = new ArrayList<>();
    list.add(getCalls(7, 8, "A"));
    list.add(getCalls(10, 11, "T"));
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
    f.filter(list, 10);
    assertFalse(list.get(0).isComplexEquivalent());
    assertFalse(list.get(1).isComplexEquivalent());
  }


}
