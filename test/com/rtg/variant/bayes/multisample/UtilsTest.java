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
