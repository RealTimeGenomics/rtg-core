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

package com.rtg.tabix;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.Resources;

import net.sf.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class TabixHeaderTest extends TestCase {

  public void test() throws IOException {
    TabixHeader th1;
    try (BlockCompressedInputStream is = new BlockCompressedInputStream(Resources.getResourceAsStream("com/rtg/tabix/resources/tabixmerge1.sam.gz.tbi"))) {
      th1 = TabixHeader.readHeader(is);
      assertEquals(4, th1.getNumSequences());
      checkOptions(th1.getOptions());
      assertTrue(Arrays.equals(new String[]{"simulatedSequence1", "simulatedSequence2", "simulatedSequence3", "simulatedSequence4"}, th1.getSequenceNamesUnpacked()));
    }
    TabixHeader th2;
    try (BlockCompressedInputStream is2 = new BlockCompressedInputStream(Resources.getResourceAsStream("com/rtg/tabix/resources/tabixmerge2.sam.gz.tbi"))) {
      th2 = TabixHeader.readHeader(is2);
      assertEquals(5, th2.getNumSequences());
      checkOptions(th2.getOptions());
      assertTrue(Arrays.equals(new String[]{"simulatedSequence4", "simulatedSequence5", "simulatedSequence6", "simulatedSequence7", "simulatedSequence8"}, th2.getSequenceNamesUnpacked()));
      final TabixHeader merged = TabixHeader.mergeHeaders(th1, th2);
      assertEquals(8, merged.getNumSequences());
      checkOptions(th2.getOptions());
      assertTrue(Arrays.equals(new String[]{"simulatedSequence1", "simulatedSequence2", "simulatedSequence3", "simulatedSequence4", "simulatedSequence5", "simulatedSequence6", "simulatedSequence7", "simulatedSequence8"}, merged.getSequenceNamesUnpacked()));
    }
  }

  private void checkOptions(TabixIndexer.TabixOptions options) {
    assertEquals(0, options.mSkip);
    assertEquals(false, options.mZeroBased);
    assertEquals('@', options.mMeta);
    assertEquals(1, options.mFormat);
    assertEquals(2, options.mSeqCol);
    assertEquals(3, options.mStartCol);
    assertEquals(-1, options.mEndCol);

  }
}
