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
package com.rtg.variant.bayes.complex;

import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;

import junit.framework.TestCase;
/**
 */
public class MendelianDenovoCheckerTest extends TestCase {
  static boolean isDenovo(DenovoChecker checker, String ref, String[] ... samples) {
    return LineageDenovoCheckerTest.isDenovo(checker, ref, samples);
  }

  public void test() {
    final Family fam = new Family("father", "mother", "child");
    final ChildFamilyLookup familyLookup = new ChildFamilyLookup(3, fam);
    final MendelianDenovoChecker checker = new MendelianDenovoChecker(familyLookup);

    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"T", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "T"}, new String[] {"T", "T"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"T", "T"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"T", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"A", "T"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"T", "T"}, new String[] {"A", "T"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"T", "T"}, new String[] {"T", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"T", "T"}, new String[] {"T", "T"}));

    assertFalse(isDenovo(checker, "T", new String[] {"A"}, new String[] {}, new String[] {"A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"A"}, new String[] {"A", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"A"}, new String[] {"T", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"A"}, new String[] {"A", "T"}, new String[] {"T", "A"}));
    assertFalse(isDenovo(checker, "T", new String[] {"T"}, null, new String[] {"T"}));
    assertFalse(isDenovo(checker, "T", null, new String[] {"T"}, new String[] {"T"}));

    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"A", "A"}));
    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"A", "C"}));
    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"}, new String[] {"T", "A"}, new String[] {"T", "C"}));

    assertTrue(isDenovo(checker, "T", new String[] {"A", "T"}, new String[] {"A"}, new String[] {"T", "T"}));
    assertTrue(isDenovo(checker, "T", new String[] {"A"}, new String[] {"A", "T"}, new String[] {"T", "T"}));
    assertTrue(isDenovo(checker, "T", new String[] {"A"}, new String[] {}, new String[] {"T"}));
    assertTrue(isDenovo(checker, "T", new String[] {"A"}, null, new String[] {"T"}));
    assertTrue(isDenovo(checker, "T", null, new String[] {"A"}, new String[] {"T"}));

    assertTrue(isDenovo(checker, "T", new String[] {"A", "A"} , new String[] {"T"}, new String[] {"T"}));
    assertTrue(isDenovo(checker, "T", new String[] {"A"} , new String[] {"T", "T"}, new String[] {"A"}));
    assertTrue(isDenovo(checker, "T", new String[] {"A", "T"} , new String[] {"T", "T"}, new String[] {"A", "A"}));
    assertTrue(isDenovo(checker, "T", new String[] {"T", "T"} , new String[] {"A", "T"}, new String[] {"A", "A"}));
  }
}
