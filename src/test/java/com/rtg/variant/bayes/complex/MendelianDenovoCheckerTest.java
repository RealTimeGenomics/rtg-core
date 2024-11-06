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
