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

package com.rtg.variant.eval;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.eval.VcfEvalParams.VcfEvalParamsBuilder;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class PathTest extends TestCase {

  public void testBestPath() {
    final byte[] template = {1, 1, 1, 1};
    final List<Variant> mutations = new ArrayList<>();
    final List<Variant> calls = new ArrayList<>();
    final MockVariant mutation = new MockVariant(2, 3, new byte[] {3}, null);
    mutations.add(mutation);
    final MockVariant call = new MockVariant(2, 3, new byte[] {3}, null);
    calls.add(call);
    final Path best = PathFinder.bestPath(template, "currentName", mutations, calls);

    assertTrue(best.getCalledIncluded().get(0).variant().equals(mutation));
    assertTrue(best.getCalledExcluded().isEmpty());
    assertTrue(best.getBaselineIncluded().get(0).variant().equals(call));
    assertTrue(best.getBaselineExcluded().isEmpty());
    final Path best2 = PathFinder.bestPath(template, "currentName", calls, mutations);
    assertTrue(best.equals(best2));
    assertTrue(best2.equals(best));
    assertFalse(best.equals(null));
    assertFalse(best.equals(calls));
    assertEquals(best.hashCode(), best2.hashCode());
  }

  private void addVar(List<OrientedVariant> list, Variant v, boolean include) {
    list.add(new OrientedVariant(v, include));
  }
  private List<Variant> getVariations(List<OrientedVariant> list) {
    final List<Variant> mutations = new ArrayList<>();
    for (final OrientedVariant v : list) {
      mutations.add(v.variant());
    }
    return mutations;
  }

  private void checkLists(List<OrientedVariant> side, List<OrientedVariant> included, List<Variant> excluded) {
    int includeCount = 0;
    int excludeCount = 0;
    for (final OrientedVariant v : side) {
      if (v.isAlleleA()) {
        includeCount++;
        //System.err.println("include: " + v);
        boolean found = false;
        for (final OrientedVariant ov : included) {
          if (ov.variant().equals(v.variant())) {
            found = true;
            break;
          }
        }
        assertTrue("the variant <" + v + "> wasn't included", found);
      } else {
        excludeCount++;
        //System.err.println("exclude: " + v);
        assertTrue("the variant <" + v.variant() + "> wasn't excluded", excluded.contains(v.variant()));
      }
    }
    assertEquals(includeCount, included.size());
    assertEquals(excludeCount, excluded.size());
  }

  private void check(byte[] template, List<OrientedVariant> aSide, List<OrientedVariant> bSide) {
    final Path best = PathFinder.bestPath(template, "currentName", getVariations(aSide), getVariations(bSide));
    //System.err.println("*****************************");
    //System.err.println(best);
    //System.err.println("*****************************");
    final List<OrientedVariant> aIncluded = best.getCalledIncluded();
    final List<Variant> aExcluded = best.getCalledExcluded();
    checkLists(aSide, aIncluded, aExcluded);

    final List<OrientedVariant> bIncluded = best.getBaselineIncluded();
    final List<Variant> bExcluded = best.getBaselineExcluded();
    checkLists(bSide, bIncluded, bExcluded);
  }

  public void testBestPath2() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(aSide, new MockVariant(4, 5, new byte[] {}, null), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(bSide, new MockVariant(5, 6, new byte[] {}, null), true);
    check(template, aSide, bSide);
  }

  public void testBestPath3() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(aSide, new MockVariant(4, 5, new byte[] {}, null), false);

    addVar(bSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    check(template, aSide, bSide);
  }

  public void testBestPath4() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, null), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(bSide, new MockVariant(5, 6, new byte[] {}, null), false);
    check(template, aSide, bSide);
  }

  public void testBestPath5() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, new byte[] {4}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {4}, new byte[] {3}), true);
    check(template, aSide, bSide);
  }

  public void testBestPath6() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, new byte[] {4}), false);

    addVar(bSide, new MockVariant(2, 3, new byte[] {4}, new byte[] {2}), false);
    check(template, aSide, bSide);
  }

  public void testBestPath7() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, new byte[] {4}), false);
    addVar(aSide, new MockVariant(5, 7, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {4}, new byte[] {2}), false);
    addVar(bSide, new MockVariant(5, 7, new byte[] {1}, new byte[] {2}), true);
    check(template, aSide, bSide);
  }

  public void testBestPath8() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(aSide, new MockVariant(3, 4, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {3}, null), true);
    addVar(bSide, new MockVariant(3, 4, new byte[] {1}, new byte[] {2}), true);
    check(template, aSide, bSide);
  }

  public void testBestPath9() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {}, null), true);
    addVar(aSide, new MockVariant(3, 4, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {}, null), true);
    addVar(bSide, new MockVariant(3, 4, new byte[] {1}, new byte[] {2}), true);
    check(template, aSide, bSide);
  }


  public void testBestPath12EdgeDeletion() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    //delete at end
    //disagee
    addVar(aSide, new MockVariant(4, 8, new byte[] {}, new byte[] {1, 1, 2}), false);
    addVar(bSide, new MockVariant(7, 8, new byte[] {3}, null), false);
    check(template, aSide, bSide);

    aSide.clear();
    bSide.clear();
    //agree
    addVar(aSide, new MockVariant(4, 8, new byte[] {}, new byte[] {1, 1, 2}), true);
    addVar(bSide, new MockVariant(4, 8, new byte[] {1, 1, 2}, new byte[] {}), true);
    check(template, aSide, bSide);

    aSide.clear();
    bSide.clear();
    //nothing in baseline
    addVar(aSide, new MockVariant(4, 8, new byte[] {}, new byte[] {1, 1, 2}), false);
    check(template, aSide, bSide);

    aSide.clear();
    bSide.clear();
    //nothing in calls
    addVar(bSide, new MockVariant(4, 8, new byte[] {1, 1, 2}, new byte[] {}), false);
    check(template, aSide, bSide);

    //disagree, both dels at edge
    aSide.clear();
    bSide.clear();
    addVar(aSide, new MockVariant(4, 8, new byte[] {}, new byte[] {1, 1, 2}), false);
    addVar(bSide, new MockVariant(6, 8, new byte[] {}, null), false);
    check(template, aSide, bSide);
  }

  public void testBestPath13EdgeInsertion() {
    //TODO currently all these have include=true due to inability to handle variants just off the edge of template
    //when fixed only the one labelled "agree" should have include=true
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    //delete at end
    //disagee
    addVar(aSide, new MockVariant(8, 8, new byte[] {}, new byte[] {3, 3, 3}), true);
    addVar(bSide, new MockVariant(8, 8, new byte[] {3}, null), true);
    check(template, aSide, bSide);
    aSide.clear();
    bSide.clear();

    //agree
    addVar(aSide, new MockVariant(8, 8, new byte[] {}, new byte[] {3, 3, 3}), true);
    addVar(bSide, new MockVariant(8, 8, new byte[] {3, 3, 3}, new byte[] {}), true);
    check(template, aSide, bSide);
    aSide.clear();
    bSide.clear();

    //only in baseline
    addVar(bSide, new MockVariant(8, 8, new byte[] {3, 3, 3}, new byte[] {}), true);
    check(template, aSide, bSide);
    aSide.clear();
    bSide.clear();

    //only in calls
    addVar(aSide, new MockVariant(8, 8, new byte[] {}, new byte[] {3, 3, 3}), true);
    check(template, aSide, bSide);
  }

  public void testBestPath10() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {3}, new byte[] {1}), false);
    addVar(aSide, new MockVariant(3, 4, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {3}, null), false);
    addVar(bSide, new MockVariant(3, 4, new byte[] {1}, new byte[] {2}), true);
    check(template, aSide, bSide);
  }
  public void testInsertion() {
    final byte[] template = {1, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 2, new byte[] {3}, null), true);

    addVar(bSide, new MockVariant(2, 2, new byte[] {3}, null), true);
    check(template, aSide, bSide);
  }

  public void testRealWorldTricky() {
    final byte[] template = {1, 1, 3, 1, 1, 1, 1, 1, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(3, 3, new byte[] {3, 1, 2, 4, 1}, new byte[] {}), true);

    addVar(bSide, new MockVariant(5, 5, new byte[] {2, 4, 1, 3, 1}, new byte[] {}), true);
    check(template, aSide, bSide);
  }
  public void testEndInclude() {
    final byte[] template = {1, 1};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(2, 3, new byte[] {2}, new byte[] {1}), true);
    check(template, aSide, bSide);
  }
  public void testEndExclude() {
    final byte[] template = {1, 1};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(2, 3, new byte[] {2}, new byte[] {1}), false);

    addVar(bSide, new MockVariant(2, 3, new byte[] {2}, new byte[] {3}), false);
    check(template, aSide, bSide);
  }

  public void testStartInclude() {
    final byte[] template = {1, 1};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 2, new byte[] {2}, new byte[] {1}), true);

    addVar(bSide, new MockVariant(1, 2, new byte[] {2}, new byte[] {1}), true);
    check(template, aSide, bSide);
  }
  public void testStartExclude() {
    final byte[] template = {1, 1};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 2, new byte[] {2}, new byte[] {1}), false);

    addVar(bSide, new MockVariant(1, 2, new byte[] {2}, new byte[] {3}), false);
    check(template, aSide, bSide);
  }

  public void testStartDelete() {
    final byte[] template = {1, 1, 3};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 2, new byte[] {}, null), true);

    addVar(bSide, new MockVariant(1, 2, new byte[] {}, null), true);
    check(template, aSide, bSide);
  }

  public void testStartDeleteExclude() {
    final byte[] template = {1, 1, 3};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 2, new byte[] {}, null), false);

    addVar(bSide, new MockVariant(1, 2, new byte[] {2}, new byte[] {}), false);
    check(template, aSide, bSide);
  }

  public void testStartInsert() {
    final byte[] template = {1, 1, 3};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 1, new byte[] {2}, null), true);

    addVar(bSide, new MockVariant(1, 1, new byte[] {2}, null), true);
    check(template, aSide, bSide);
  }

  public void testStartInsertExclude() {
    final byte[] template = {1, 1, 3};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, new MockVariant(1, 1, new byte[] {2}, null), false);

    addVar(bSide, new MockVariant(1, 1, new byte[] {3}, null), false);
    check(template, aSide, bSide);
  }

//  public void testOverlappingError() {
//    Diagnostic.setLogStream(SimpleTestUtils.getNullPrintStream());
//    try {
//      byte[] template = {1, 1, 3, 1, 2};
//      List<OrientedVariant> aSide = new ArrayList<>();
//      List<OrientedVariant> bSide = new ArrayList<>();
//      addVar(aSide, new MockVariant(1, 3, new byte[] {2}, null), false);
//      addVar(aSide, new MockVariant(2, 4, new byte[] {2}, null), false);
//
//      addVar(bSide, new MockVariant(1, 1, new byte[] {3}, null), false);
//      try {
//        check(template, aSide, bSide);
//        fail("Expected an exception");
//      } catch (NoTalkbackSlimException e) {
//        // Expected
//        assertTrue(e.getMessage().contains("Overlapping variants aren't supported"));
//      }
//    } finally {
//      Diagnostic.setLogStream();
//    }
//  }

  private DetectedVariant getDetectedVariant(String vcfLine) {
    final String linetab = vcfLine.replaceAll(" ", "\t");
    return new DetectedVariant(VcfReader.vcfLineToRecord(linetab), 0, RocSortValueExtractor.NULL_EXTRACTOR, false);
  }

  public void testWithDetectedVariant() {
    final String line = "simulatedSequence1 3 . C T,G 0.0 PASS . GT 1/2";
    final DetectedVariant v = getDetectedVariant(line);
    final String line2 = "simulatedSequence1 3 . C G,T 31.0 PASS . GT 1/2";
    final DetectedVariant v2 = getDetectedVariant(line2);

    final byte[] template = {1, 1, 2, 4, 4};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, v2, true);

    addVar(bSide, v, true);
    check(template, aSide, bSide);
  }

  public void testDontCrossTheStreams() {
    final DetectedVariant vA = getDetectedVariant("simulatedSequence1 1 . A T,G 0.0 PASS . GT 1/2");
    final DetectedVariant vA2 = getDetectedVariant("simulatedSequence1 1 . A G,T 0.0 PASS . GT 1/2");
    final DetectedVariant v = getDetectedVariant("simulatedSequence1 3 . C T,G 0.0 PASS . GT 1/2");
    final DetectedVariant v2 = getDetectedVariant("simulatedSequence1 3 . C G,T 31.0 PASS . GT 1/2");

    final byte[] template = {1, 1, 2, 4, 4};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, vA, true);
    addVar(aSide, v2, true);

    addVar(bSide, vA2, true);
    addVar(bSide, v, true);
    check(template, aSide, bSide);
  }

  public void testRealWordTricky() {
//    final String line = "chr22 5 e ga.a CTAGA:i 25.5 42 4.578 GAA 19 2.071 GACTAGAA 18 1.962 GTCTAGAA 1 0.109 G~ 1 0.109 TACTAGAA 1 0.109 ~A 2 0.218".replaceAll(" ", "\t");
    final String line = "chr22 4 . A ACTAGA 25.5 PASS . GT 0/1";
//    final String line2 = "chr22 3 e i GACTA:i".replaceAll(" ", "\t");
    final String line2 = "chr22 2 . A AGACTA 0.0 PASS . GT 0/1";
    final DetectedVariant v1 = getDetectedVariant(line);
    final DetectedVariant v2 = getDetectedVariant(line2);

    final byte[] template = {1, 1, 3, 1, 1, 2, 2};
    final List<OrientedVariant> aSide = new ArrayList<>();
    final List<OrientedVariant> bSide = new ArrayList<>();
    addVar(aSide, v1, true);
    addVar(bSide, v2, true);
    check(template, aSide, bSide);
  }

  public void testHeterozygousInsert()  {
    final byte[] template = {1, 2, 1, 4, 4, 4, 4, 1, 2,  1};
    final OrientedVariant[] a = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 7 . T TT 0.0 PASS . GT 0/1"), false)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final OrientedVariant[] b = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 4 . T TT 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final double[] expectedWeights = {1.0, 1.0};

    Path original = PathFinder.bestPath(template, "currentName", Arrays.asList(a), Arrays.asList(b));
    Path.calculateWeights(original, original.getCalledIncluded(), original.getBaselineIncluded());
    check(original.getCalledIncluded(), expectedWeights);

    Path originalRev = PathFinder.bestPath(template, "currentName", Arrays.asList(b), Arrays.asList(a));
    Path.calculateWeights(originalRev, originalRev.getCalledIncluded(), originalRev.getBaselineIncluded());
    check(originalRev.getCalledIncluded(), expectedWeights);
  }

  static void check(List<OrientedVariant> v, double[] weights) {
    assertEquals(v.size(), weights.length);
    for (int i = 0; i < v.size(); i++) {
      assertEquals(i + " : " + v.toString(), weights[i], v.get(i).getWeight(), 0.0001);
    }
  }

  public void testHeterozygousDelete()  {
    final byte[] template = {1, 2, 1, 4, 4, 4, 4, 1, 2, 1};
    final OrientedVariant[] a = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 6 . TT T 0.0 PASS . GT 1/0"), false)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final OrientedVariant[] b = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 4 . TT T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final double[] expectedWeights = {1.0, 1.0};

    Path original = PathFinder.bestPath(template, "currentName", Arrays.asList(a), Arrays.asList(b));
    Path.calculateWeights(original, original.getCalledIncluded(), original.getBaselineIncluded());
    check(original.getCalledIncluded(), expectedWeights);

    Path originalRev = PathFinder.bestPath(template, "currentName", Arrays.asList(b), Arrays.asList(a));
    Path.calculateWeights(originalRev, originalRev.getCalledIncluded(), originalRev.getBaselineIncluded());
    check(originalRev.getCalledIncluded(), expectedWeights);
  }

  public void testNop()  {
    final byte[] template = {1, 2, 1, 4, 4, 4, 4, 1, 2, 1};
    final OrientedVariant[] a = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 3 . AT A 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 6 . T TT 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final OrientedVariant[] b = {
        new OrientedVariant(getDetectedVariant("seq 2 . C T 0.0 PASS . GT 1/1"), true)
      , new OrientedVariant(getDetectedVariant("seq 9 . C T 0.0 PASS . GT 1/1"), true)
    };

    final double[] expectedWeights = {0.3333, 0.3333, 0.3333, 1.0};

    Path original = PathFinder.bestPath(template, "currentName", Arrays.asList(a), Arrays.asList(b));
    Path.calculateWeights(original, original.getCalledIncluded(), original.getBaselineIncluded());
    check(original.getCalledIncluded(), expectedWeights);
  }
  
  public void testPastEndOfTemplate() {
    final byte[] template = {1, 2, 1, 4, 4, 4, 4, 1, 2, 1};
    final OrientedVariant[] a = {new OrientedVariant(getDetectedVariant("seq 12 . C T 0.0 PASS . GT 1/1"), true)};
    final List<OrientedVariant> al = Arrays.asList(a);
    //Testing that having variants outside the template will not attempt to advance the path past the end of the
    //template when the path is in sync
    PathFinder.bestPath(template, "seq", al, al);
  }

  public void testNoInfiniloop() throws Exception {
    final String templatefa = ">CFTR.8.500s\nTCATCACTAAGGTTAGCATGTAATAGTACAAGGAAGAATCAGTTGTATGTTAAATCTAATGTATAAAAAGTTTTATAAAATATCATATGTTTAGAGAGTATATTTCAAATATGTTGAATCCTAGTGCTTGGGTGCAAATTAACTTTAGAACACTAGTAAAATTATTTTATTAAGAAATAATTACTATTTCATTATTAAAATTCATATATAAGATGTAGCACAATGAGAGTATAAAGTAGATGTAATAATGCATTAATGCTATTCTGATTCTATAATATGTTTTTGCTCTCTTTTATAAATAGGATTTCTTACAAAAGCAAGAATATAAGACATTGGAATATAACTTAACGACTACAGAAGTAGTGATGGAGAATGTAACAGCCTTCTGGGAGGAGGTCAGATAATTTTTAAAAAATTGTTTGCTCTAAACACCTAACTGTTTTCTTCTTTGTGAATATGGCCTAATGGCGAATAAAATTAGAATGATGATATAACTGGTAGAACTGGAAGGAGGATCACT\n";

    final String generated = VcfHeader.MINIMAL_HEADER + "\tSAMPLE\n"
        + "CFTR.8.500s\t114\t.\tT\tA\t.\tPASS\t.\tGT:GQ\t1/1:3.0\n";

    final String called = VcfHeader.MINIMAL_HEADER + "\tCD0000\n"
        + "CFTR.8.500s\t114\t.\tT\tA\t222\t.\tDP=510;VDB=0.0090;AF1=1;AC1=2;DP4=3,0,287,220;MQ=53;FQ=-282;PV4=0.26,1,0.37,1\tGT:PL:GQ\t1/1:255,255,0:99\n"
        + "CFTR.8.500s\t518\t.\tACT\tACTCACTTATTTTCT\t11.8\t.\tINDEL;DP=185;VDB=0.0114;AF1=1;AC1=2;DP4=130,40,0,15;MQ=43;FQ=-290;PV4=2.7e-09,1,1.6e-05,0.0017\tGT:PL:GQ\t1/1:52,255,0:75\n";

    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File template = ReaderTestUtils.getDNADir(templatefa, new File(tmpDir, "template"));
      final File muts = FileHelper.stringToGzFile(generated, new File(tmpDir, "gen.vcf.gz"));
      new TabixIndexer(muts).saveVcfIndex();
      final File calls = FileHelper.stringToGzFile(called, new File(tmpDir, "calls.vcf.gz"));
      new TabixIndexer(calls).saveVcfIndex();

      final File output = new File(tmpDir, "eval");

      final OutputParams op = new OutputParams(output, false, false);
      final VcfEvalParams mep = new VcfEvalParamsBuilder().templateFile(template).baseLineFile(muts).callsFile(calls).outputParams(op).useAllRecords(true).create();
      VcfEvalTask.evaluateCalls(mep);

      final String weightedroc = FileUtils.fileToString(new File(output, "weighted_roc.tsv"));
      TestUtils.containsAll(weightedroc, "99.000\t1.000\t0" + StringUtils.LS, "75.000\t1.000\t1" + StringUtils.LS);
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
