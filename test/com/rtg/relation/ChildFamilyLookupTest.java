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
package com.rtg.relation;

import junit.framework.TestCase;

/**
 */
public class ChildFamilyLookupTest extends TestCase {

  public void test() {
    final Family f1 = new Family("father1", "mother", "child1");
    final Family f2 = new Family("father2", "mother", "child2");
    f2.setSampleId(1, 1);
    f2.setSampleId(0, 3);
    f2.setSampleId(2, 4);
    final Family[] fams = {f1, f2};
    assertEquals(2, fams.length);
    final ChildFamilyLookup lookup = new ChildFamilyLookup(5, fams);
    for (final Family f : fams) {
      final int[] sampleIds = f.getSampleIds();
      assertNull(lookup.getFamily(sampleIds[Family.FATHER_INDEX]));
      assertNull(lookup.getFamily(sampleIds[Family.MOTHER_INDEX]));
      for (int i = Family.FIRST_CHILD_INDEX; i < sampleIds.length; i++) {
        assertEquals(f, lookup.getFamily(sampleIds[i]));
      }
    }
  }
}
