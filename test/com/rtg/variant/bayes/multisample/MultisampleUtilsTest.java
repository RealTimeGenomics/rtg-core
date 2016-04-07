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

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Covariate;
import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.CalibratedMachineErrorChooser;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.ReadGroupMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.HypothesesComplex;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class MultisampleUtilsTest extends TestCase {

  public void testUserSelected() throws Exception {
    Diagnostic.setLogStream();
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.machineErrorName("illumina");
    final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(b.create());
    assertTrue(chooser instanceof DefaultMachineErrorChooser);
  }

  public void testUncalibrated() throws Exception {
    Diagnostic.setLogStream();
    final ArrayList<File> mapped = new ArrayList<>();
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.mapped(mapped);
    final SAMFileHeader header = new SAMFileHeader();
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
    rg.setPlatform("illumina");
    header.addReadGroup(rg);
    b.uberHeader(header);
    final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(b.create());
    assertTrue(chooser instanceof ReadGroupMachineErrorChooser);
  }

  public void testCalibrated() throws Exception {
    Diagnostic.setLogStream();
    final ArrayList<File> mapped = new ArrayList<>();
    final Calibrator c = new Calibrator(new Covariate[0], null);
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.mapped(mapped);
    b.calibrator(c);
    final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(b.create());
    assertTrue(chooser instanceof CalibratedMachineErrorChooser);
  }

  private static VariantAlignmentRecord makeAlignment(final int start, final String read, final String cigar) {
    return makeAlignment(start, start, read, cigar);
  }

  private static VariantAlignmentRecord makeAlignment(final int id, final int start, final String read, final String cigar) {
    final SAMRecord s1 = makeSamRecord(start, read, cigar);
    s1.setReadName("r" + id);
    s1.setMappingQuality(30);
    return new VariantAlignmentRecord(s1);
  }

  public static SAMRecord makeSamRecord(final int alignStart, final String readString, final String cigar) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(alignStart);
    sam.setReadString(readString);
    sam.setCigarString(cigar);
    sam.setBaseQualities(new byte[readString.length()]);
    return sam;
  }
  /*
   *    1234567_890123456
   *    acgtgtg_gtcgtacgtaccaagtaca
   *    acgtgtg
   *     cgtgtga
   *      gtgtgcg
   *           ggtcgta
   *            gtcgtac
   * Simple(!) insert
   */
  public void testIntersectSet1() throws Exception {
    final TreeSet<VariantAlignmentRecord> set = new TreeSet<>();
    set.add(makeAlignment(1, "ACGTGTG", "7="));
    set.add(makeAlignment(2, "CGTGTGA", "6=1I"));
    set.add(makeAlignment(3, "GTGTGCG", "5=1I1="));
    set.add(makeAlignment(8, "GGTCGTA", "1I6="));
    set.add(makeAlignment(8, "GTCGTAC", "7="));
    final List<AlignmentMatch> is = MultisampleUtils.intersectSet(8 - 1, 8 - 1, set.iterator(), null, new VariantParamsBuilder().create());
    //System.err.println(is);
    assertEquals(3, is.size());
    assertEquals("A", is.get(0).toString());
    assertEquals("C", is.get(1).toString());
    assertEquals("G", is.get(2).toString());

    final VariantParams p = VariantParams.builder().genomePriors("human").create();
    ComplexTemplate context = new ComplexTemplate(DnaUtils.encodeString("acgtgtggtcgtacgtaccaagtaca"), "", 8 - 1, 8 - 1);
    context.setComplexContext(HypothesesComplex.createComplexDescription(is, context, null, p.pruneHypotheses(), p.maxComplexHypotheses()), SimplePossibility.SINGLETON);
    final HypothesesComplex hc = HypothesesComplex.makeComplexHypotheses(context, false, p);
    final HashSet<String> haphypoexp = new HashSet<>();
    haphypoexp.add("");
    haphypoexp.add("A");
    haphypoexp.add("C");
    haphypoexp.add("G");
    int hap = 0;
    int het = 0;
    final HashSet<String> haphypoact = new HashSet<>();
    for (int i = 0; i < hc.size(); i++) {
      if (hc.code().homozygous(i)) {
        hap++;
        haphypoact.add(hc.description().name(i));
      } else {
        het++;
      }
    }
    assertEquals(haphypoexp, haphypoact);
    assertEquals(4, hap); //i, A, C, G
    assertEquals(6, het); //i:A, i:C, i:G, A:C, A:G, C:G
    assertEquals("", hc.description().name(hc.reference()));
  }

  /*
   *    1234567_890_123456
   *    acgtgtg_gtc_gtacgtaccaagtaca
   *     cgtgtgc
   *      gtgtg_gt
   *             tc_gtacg
   *               agtacgt
   */
  public void testIntersectSet2() throws Exception {
    final TreeSet<VariantAlignmentRecord> set = new TreeSet<>();
    set.add(makeAlignment(2, "CGTGTGC", "6=1I"));
    set.add(makeAlignment(3, "GTGTGGT", "7="));
    set.add(makeAlignment(9, "TCGTACG", "7="));
    set.add(makeAlignment(11, "AGTACGT", "1I6="));
    final List<AlignmentMatch> is = MultisampleUtils.intersectSet(8 - 1, 11 - 1, set.iterator(), null, new VariantParamsBuilder().create());
    //System.err.println(is);
    assertEquals(4, is.size());
    assertEquals("C~", is.get(0).toString());
    assertEquals("GT~", is.get(1).toString());
    assertEquals("~TC", is.get(2).toString());
    assertEquals("~A", is.get(3).toString());
    final VariantParams p = VariantParams.builder().genomePriors("human").create();
    ComplexTemplate context = new ComplexTemplate(DnaUtils.encodeString("acgtgtggtcgtacgtaccaagtaca"), "", 8 - 1, 11 - 1);
    context.setComplexContext(HypothesesComplex.createComplexDescription(is, context, null, p.pruneHypotheses(), p.maxComplexHypotheses()), SimplePossibility.SINGLETON);
    final HypothesesComplex hc = HypothesesComplex.makeComplexHypotheses(context, false, p);
    assertEquals(1, hc.size()); //just the ref, rest are unfixed
    assertEquals("GTC", hc.description().name(hc.reference()));
  }
  /*
   *    1234567890123456
   *    acgtgtggtcgtacgtaccaagtaca
   *     cgtgtgg
   *     cgtg*gg
   *     cgtg*gg
   *     cgtg*gg
   * simple delete
   */
  public void testIntersectSet3() throws Exception {
    final TreeSet<VariantAlignmentRecord> set = new TreeSet<>();
    set.add(makeAlignment(1, 2, "CGTGTGG", "7="));
    set.add(makeAlignment(2, 2, "GGTGGG", "1X3=1D2="));
    set.add(makeAlignment(3, 2, "CGTGGG", "4=1D2="));
    set.add(makeAlignment(4, 2, "CGTGGG", "4=1D2="));
    final List<AlignmentMatch> is = MultisampleUtils.intersectSet(6 - 1, 7 - 1, set.iterator(), null, new VariantParamsBuilder().create());
    assertEquals(4, is.size());
    assertEquals("T", is.get(0).toString());
    assertEquals("", is.get(1).toString());
    assertEquals("", is.get(2).toString());
    assertEquals("", is.get(3).toString());
    final VariantParams p = VariantParams.builder().genomePriors("human").create();
    ComplexTemplate context = new ComplexTemplate(DnaUtils.encodeString("acgtgtggtcgtacgtaccaagtaca"), "", 6 - 1, 7 - 1);
    context.setComplexContext(HypothesesComplex.createComplexDescription(is, context, null, p.pruneHypotheses(), p.maxComplexHypotheses()), SimplePossibility.SINGLETON);
    final HypothesesComplex hc = HypothesesComplex.makeComplexHypotheses(context, false, p);
    assertEquals(3, hc.size());
    final HashSet<String> haphypoexp = new HashSet<>();
    haphypoexp.add("");
    haphypoexp.add("T");
    int hap = 0;
    int het = 0;
    final HashSet<String> haphypoact = new HashSet<>();
    for (int i = 0; i < hc.size(); i++) {
      if (hc.code().homozygous(i)) {
        hap++;
        haphypoact.add(hc.description().name(i));
      } else {
        het++;
      }
    }
    assertEquals(haphypoexp, haphypoact);
    assertEquals(2, hap);
    assertEquals(1, het);
    assertEquals("T", hc.description().name(hc.reference()));
  }


}
