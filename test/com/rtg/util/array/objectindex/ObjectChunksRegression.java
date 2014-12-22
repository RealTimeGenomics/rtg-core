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

package com.rtg.util.array.objectindex;

import com.rtg.util.test.RandomByteGenerator;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Class for testing the Object Chunks implementation that is capable of
 * holding more than the maximum integer worth of data.
 */
public class ObjectChunksRegression extends TestCase {

  private static final long NUM_ELEMENTS = 2L * Integer.MAX_VALUE + 9000L;
  private static final int RANGE = 128;

  /**
   * Test the common index implementation
   */
  public void testIndex() {
    Diagnostic.setLogStream();
    doTest(RANGE, NUM_ELEMENTS);
  }

  private void doTest(int range, long elements) {
    final RandomByteGenerator value = new RandomByteGenerator(range);

    final ObjectChunks<Byte> index = new ObjectChunks<>(elements);
    assertEquals(elements, index.length());

    for (long l = 0; l < elements; l++) {
      index.set(l, value.nextValue());
    }

    value.reset();

    for (long l = 0; l < elements; l++) {
      assertEquals(Byte.valueOf(value.nextValue()), index.get(l));
    }
  }
}
