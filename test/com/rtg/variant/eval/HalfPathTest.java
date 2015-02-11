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

import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.BasicLinkedListNode;
import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfReader;

import junit.framework.TestCase;

/**
 * Test the corresponding class
 */
public class HalfPathTest extends TestCase {

  private static <T> List<T> asList(BasicLinkedListNode<T> vals) {
    final ArrayList<T> list = new ArrayList<>();
    for (T val : vals) {
      list.add(val);
    }
    return list;
  }


  public void testNextBase() {
    //                            1  2  3  4  5  6  7  8  9
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1};
                           //     2, 1  1 31  1     2  1  1
    final HalfPath path = new HalfPath(template);
    assertEquals(0, path.getVariantEndPosition());
    final OrientedVariant included = new OrientedVariant(new MockVariant(1, 2, new byte[] {2}, null), true);
    path.include(included);
    assertEquals(1, path.getVariantEndPosition());
    path.include(new OrientedVariant(new MockVariant(4, 4, new byte[] {3}, null), true));
    path.include(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, null), true));

    final MockVariant excluded = new MockVariant(3, 4, new byte[] {3}, null);
    path.exclude(excluded);
    assertEquals(3, path.getVariantEndPosition());
    assertTrue(asList(path.getIncluded()).contains(included));
    assertTrue(asList(path.getExcluded()).contains(excluded));

    final byte[] expected = {2, 1, 1, 3, 1, 1, 2, 1, 1};
    int i = 0;
    path.step();
    while (!path.finished()) {
      assertEquals(expected[i], path.nextHaplotypeABase());
      path.step();
      i++;
    }

  }

  public void testHetero() {
    //                            1  2  3  4  5  6  7  8  9
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1};
                           //     2, 1  1 31  1     2  1  1
                           //     2, 1  1431  1     2  1  1
    final HalfPath path = new HalfPath(template);
    path.include(new OrientedVariant(new MockVariant(1, 2, new byte[] {2}, new byte[] {1}), false));
    path.include(new OrientedVariant(new MockVariant(4, 4, new byte[] {3}, new byte[] {4, 3}), true));
    path.include(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, null), true));

    final byte[] expected = {1, 1, 1, 3, 1, 1, 2, 1, 1};
    final byte[] expectedMinus = {2, 1, 1, 4, 3, 1, 1, 2, 1, 1};
    int i = 0;
    path.step();
    while (!path.finished()) {
      if (!path.finishedHaplotypeA()) {
        assertEquals(expected[i], path.nextHaplotypeABase());
      }
      if (!path.finishedHaplotypeB()) {
        assertEquals(expectedMinus[i], path.nextHaplotypeBBase());
        //System.out.println(expectedMinus[i] + " " + path.nextMinusBase());
      }
      path.step();
      i++;
    }
    assertEquals(Math.max(expected.length, expectedMinus.length), i);

  }

  public void testCompare() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    assertFalse(path.equals(null));
    assertFalse(path.finishedHaplotypeA());
    assertFalse(path.finishedHaplotypeB());
    assertEquals(0, path.getVariantEndPosition());

    path.include(new OrientedVariant(new MockVariant(1, 2, new byte[] {2}, new byte[] {1}), false));
    assertEquals(1, path.getVariantEndPosition());
    path.include(new OrientedVariant(new MockVariant(4, 4, new byte[] {3}, new byte[] {4, 3}), true));
    path.include(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, null), true));
    final HalfPath copy = new HalfPath(path);
    assertEquals(0, path.compareTo(copy));
    final OrientedVariant var = new OrientedVariant(new MockVariant(10, 11, new byte[] {2}, null), true);
    final OrientedVariant reverse = new OrientedVariant(new MockVariant(10, 11, new byte[] {3}, new byte[] {2}), false);

    copy.include(var);
    assertEquals(-1, path.compareTo(copy));
    assertEquals(1, copy.compareTo(path));

    path.step();
    assertEquals(1, path.compareTo(copy));
    assertEquals(-1, copy.compareTo(path));

    copy.step();
    assertEquals(-1, path.compareTo(copy));
    assertEquals(1, copy.compareTo(path));

    path.include(reverse);
    assertEquals(-1, path.compareTo(copy));
    assertEquals(1, copy.compareTo(path));

    while (path.getPosition() < 8) {
      path.step();
    }
    while (copy.getPosition() < 8) {
      copy.step();
    }
    assertEquals(-1, path.compareTo(copy));
    path.step();
    copy.step();
    assertEquals(-1, path.compareTo(copy));
    final OrientedVariant both = new OrientedVariant(new MockVariant(12, 12, new byte[] {1, 2, 3}, null), true);
    path.include(both);
    copy.include(both);
    while (path.getPosition() < 12) {
      path.step();
    }
    while (copy.getPosition() < 12) {
      copy.step();
    }
    assertEquals(0, copy.compareTo(path));
    assertEquals(0, path.compareTo(copy));
    assertTrue(copy.hashCode() == path.hashCode());
    assertTrue(copy.equals(path));
    copy.step();
    assertFalse(copy.hashCode() == path.hashCode());
    assertFalse(copy.equals(path));
    assertEquals(-1, path.compareTo(copy));
    assertEquals(1, copy.compareTo(path));
  }

  public void testMoreCompares() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    final HalfPath copy = new HalfPath(template);
    path.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, new byte[] {4, 3}), true));
    copy.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, new byte[] {4, 3}), true));
    path.step();
    path.step();
    copy.step();
    assertEquals(-1, copy.compareTo(path));
    assertEquals(1, path.compareTo(copy));
  }

  public void testMoreCompares2() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    final HalfPath copy = new HalfPath(template);
    path.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, null), true));
    copy.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, new byte[] {4, 3}), true));
    assertEquals(1, copy.compareTo(path));
    assertEquals(-1, path.compareTo(copy));
  }

  public void testMoreCompares3() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    final HalfPath copy = new HalfPath(template);
    path.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, null), true));
    copy.include(new OrientedVariant(new MockVariant(2, 3, new byte[] {3, 4}, null), true));
    assertEquals(0, copy.compareTo(path));
  }

  public void testPosition() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    path.step();
    path.step();
    assertEquals(1, path.getPosition());
    path.include(new OrientedVariant(new MockVariant(3, 4, new byte[] {2, 3, 4}, new byte[] {4}), true));
    path.step();
    path.step();
    path.step();
    assertEquals(4, path.getPosition());

    final HalfPath path2 = new HalfPath(template);
    path2.step();
    path2.step();
    assertEquals(1, path2.getPosition());
    path2.include(new OrientedVariant(new MockVariant(3, 4, new byte[] {4}, new byte[] {2, 3, 4}), true));
    path2.step();
    path2.step();
    path2.step();
    assertEquals(4, path2.getPosition());
  }

  public void testToString() {
    final byte[] template = {1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4};
    final HalfPath path = new HalfPath(template);
    assertFalse(path.toString().contains("minus:"));
    path.include(new OrientedVariant(new MockVariant(3, 4, new byte[] {3}, new byte[] {4}), true));
    path.exclude(new MockVariant(5, 6, new byte[] {2}, new byte[] {4}));
    path.exclude(new MockVariant(5, 6, new byte[] {2}, new byte[] {4}));
    path.include(new OrientedVariant(new MockVariant(8, 9, new byte[] {3}, new byte[] {4}), true));
    TestUtils.containsAll(path.toString()
        , "included:"
        , "[8:9 G:T+, 3:4 G:T+]"
        , "excluded:"
        , "[5:6 C:T, 5:6 C:T]");

  }
  private static final String[] CALLS_TRICKY = {
      "seq 1 . A T 1 PASS . GT 1/1"
    , "seq 2 . A AT 9 PASS . GT 1/1"
    , "seq 4 . GT G 2 PASS . GT 1/1"
    , "seq 10 . CGT AGA 3 PASS . GT 1/1"
    , "seq 14 . C A 8 PASS . GT 1/1"
    , "seq 16 . T A 4 PASS . GT 1/1"
    , "seq 20 . C A 10 PASS . GT 1/1 "
  };
  private static final String[] MUTATIONS_TRICKY = {
      "seq 5 . TC T 0.0 PASS . GT 1/1"
    , "seq 7 . G GC 0.0 PASS . GT 1/1"
    , "seq 10 . C A 0.0 PASS . GT 1/1"
    , "seq 12 . T A 0.0 PASS . GT 1/1"
    , "seq 14 . CGT AGA 0.0 PASS . GT 1/1"
    , "seq 18 . C A 0.0 PASS . GT 1/1"
    , "seq 20 . C A 0.0 PASS . GT 1/1"

  };

  public void testCaseFromMutationEvalTricky() {
    final byte[] template = DnaUtils.encodeString("ACTTTCCCACGTACGTCCTCT");
    HalfPath path = new HalfPath(template);

    for (final String var : CALLS_TRICKY) {
      final String vartab = var.replaceAll(" ", "\t");
      path.include(new OrientedVariant(new DetectedVariant(VcfReader.vcfLineToRecord(vartab), 0, RocSortValueExtractor.NULL_EXTRACTOR, false), true));
    }
    StringBuilder sb = new StringBuilder();
    path.step();
    while (!path.finished()) {
      sb.append(DNA.valueChars()[path.nextHaplotypeABase()]);
      path.step();
    }
    assertEquals("TCTTTCCCAAGAAAGACCTAT", sb.toString());

    path = new HalfPath(template);
    for (final String var : MUTATIONS_TRICKY) {
      final String vartab = var.replaceAll(" ", "\t");
      path.include(new OrientedVariant(new DetectedVariant(VcfReader.vcfLineToRecord(vartab), 0, RocSortValueExtractor.NULL_EXTRACTOR, false), true));
    }
    sb = new StringBuilder();
    path.step();
    while (!path.finished()) {
      sb.append(DNA.valueChars()[path.nextHaplotypeABase()]);
      path.step();
    }
    assertEquals("ACTTTCCCAAGAAAGACATAT", sb.toString());
  }
}
