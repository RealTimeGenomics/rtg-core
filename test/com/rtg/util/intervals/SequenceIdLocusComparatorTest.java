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
package com.rtg.util.intervals;

import junit.framework.TestCase;

/**
 */
public class SequenceIdLocusComparatorTest extends TestCase {


  public void testCompare() {
    SequenceIdLocus pos1 = new SequenceIdLocusSimple(10, 100, 102);
    SequenceIdLocusComparator comp = new SequenceIdLocusComparator();

    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(11, 99, 101)) < 0);
    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(9, 99, 101)) > 0);

    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(10, 100, 102)) == 0);
    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(10, 101, 103)) < 0);
    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(10, 99, 101)) > 0);

    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(10, 100, 110)) < 0);
    assertTrue(comp.compare(pos1, new SequenceIdLocusSimple(10, 100, 101)) > 0);

  }

}
