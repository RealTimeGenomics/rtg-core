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

package com.rtg.reader;

import com.rtg.util.test.RandomByteGenerator;
import com.rtg.util.bytecompression.ByteArray;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Base test class for testing the various <code>ByteArray</code> implementations
 * that can have large amounts of data contained within them.
 */
public abstract class AbstractByteArrayRegression extends TestCase {

  private static final int RANGE = 128;
  private static final long NUM_ELEMENTS = 10L * Integer.MAX_VALUE + 9000L;

  protected abstract ByteArray createByteArray(int range, long elements);

  /**
   * Test the byte array
   */
  public void testByteArrayImplementation() {
    Diagnostic.setLogStream();
    doTest(RANGE, NUM_ELEMENTS);
  }

  private void doTest(int range, long elements) {
    final RandomByteGenerator value = new RandomByteGenerator(range);
    final byte[] buffer = new byte[1024];
    final ByteArray byteArray = createByteArray(range, elements);
    assertEquals(elements, byteArray.length());
    long lastOffset = 0;
    for (long l = 0; l < elements;) {
      int i = 0;
      for (; i < buffer.length && l < elements; i++, l++) {
        buffer[i] = value.nextValue();
      }
      byteArray.set(lastOffset, buffer, i);
      lastOffset = l;
    }

    value.reset();

//  final long labelChunkSize = NUM_ELEMENTS / 1000L;
//  long nextLabelOutput = labelChunkSize;
    byte val = -1;
    for (long l = 0; l < elements;) {
      final int read = (int) Math.min(buffer.length, elements - l);
      byteArray.get(buffer, l, read);
      for (int i = 0; i < read; i++, l++) {
        val = value.nextValue();
        assertEquals(val, buffer[i]);
      }
//      if (l >= nextLabelOutput) {
//        System.err.println("Elements Read: " + l);
//        nextLabelOutput += labelChunkSize;
//      }
    }
    assertEquals(val, byteArray.get(elements - 1));
  }
}
