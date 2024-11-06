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

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.relation.Family;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class UtilsTest extends TestCase {

  public void testSexMemoNull() throws IOException {
    Diagnostic.setLogStream();
    try {
      Utils.createSexMemo(VariantParams.builder().create());
      fail("expected NPE");
    } catch (NullPointerException e) {
      // Expected
    }
  }

  public void testNumberAlleles() {
    final Code code = new CodeDiploid(4);
    assertEquals(3, Utils.numberAlleles(4, code, 0, code.code(1, 2)));
    assertEquals(0, Utils.numberAlleles(4, code));
    assertEquals(1, Utils.numberAlleles(4, code, 0));
    assertEquals(2, Utils.numberAlleles(4, code, 0, 1));
    assertEquals(2, Utils.numberAlleles(4, code, 1, 0, 0, 1));
    assertEquals(4, Utils.numberAlleles(4, code, 0, code.code(1, 2), code.code(2, 3)));
    assertEquals(4, Utils.numberAlleles(4, code, code.code(1, 2), code.code(2, 3), 0));
    assertEquals(4, Utils.numberAlleles(4, code, code.code(1, 2), 0, code.code(2, 3)));
  }

  public void testNumberAllelesExcl() {
    final Code code = new CodeDiploid(4);
    assertEquals(2, Utils.numberAllelesExclude(4, code, 0, code.code(1, 2)));
    assertEquals(0, Utils.numberAllelesExclude(4, code, 0));
    assertEquals(1, Utils.numberAllelesExclude(4, code, 0, 1));
    assertEquals(2, Utils.numberAllelesExclude(4, code, 0, 1, 2));
    assertEquals(0, Utils.numberAllelesExclude(4, code, 0, 0));
    assertEquals(1, Utils.numberAllelesExclude(4, code, 1, 0, 0, 1));
    assertEquals(1, Utils.numberAllelesExclude(4, code, 1, 0, 1));
    assertEquals(0, Utils.numberAllelesExclude(4, code, 1, 1, 1));
    assertEquals(2, Utils.numberAllelesExclude(4, code, 0, code.code(0, 2), code.code(2, 3)));
    assertEquals(2, Utils.numberAllelesExclude(4, code, 0, code.code(2, 3), code.code(0, 2)));

    assertEquals(3, Utils.numberAllelesExclude(4, code, -1, 0, code.code(1, 2)));
    assertEquals(0, Utils.numberAllelesExclude(4, code, -1));
    assertEquals(1, Utils.numberAllelesExclude(4, code, -1, 0));
    assertEquals(2, Utils.numberAllelesExclude(4, code, -1, 0, 1));
    assertEquals(2, Utils.numberAllelesExclude(4, code, -1, 1, 0, 0, 1));
    assertEquals(4, Utils.numberAllelesExclude(4, code, -1, 0, code.code(1, 2), code.code(2, 3)));
    assertEquals(4, Utils.numberAllelesExclude(4, code, -1, code.code(1, 2), code.code(2, 3), 0));
    assertEquals(4, Utils.numberAllelesExclude(4, code, -1, code.code(1, 2), 0, code.code(2, 3)));

  }

  public void testIsCallableAsFamily() {
    final Family family = new Family("father", "mother", "child1", "child2", "child3");

    final ArrayList<String> outputs = new ArrayList<>();
    assertFalse(Utils.isCallableAsFamily(outputs, family));

    outputs.add("father");
    assertFalse(Utils.isCallableAsFamily(outputs, family));

    outputs.add("mother");
    assertFalse(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child3");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child2");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child1");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("other");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.clear();
    outputs.add("father");
    assertFalse(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child2");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child3");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.clear();
    outputs.add("child1");
    assertFalse(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child2");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("child3");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

    outputs.add("mother");
    assertTrue(Utils.isCallableAsFamily(outputs, family));

  }
}
