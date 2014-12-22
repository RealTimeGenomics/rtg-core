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
import java.io.InputStream;

import com.rtg.util.Resources;

import net.sf.samtools.util.BlockCompressedInputStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SamPositionReaderTest extends TestCase {

  private static final String[] EXP_REF_NAME = {"gi|89161203|ref|NC_000022.9|NC_000022", "*"};
  private static final int[] ENTRIES = {15, 15};
  private static final int[][] START = {new int[] {14430018, 14430080, 14430108, 14430124, 14430131,
                                        14430133, 14430187, 14430210, 14430243, 14430246,
                                        14430251, 14430255, 14430302, 14430309, 14430328},
                                        new int[] {0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0}
  };

  private static final int[][] LENGTH = {new int[] {100, 100, 100, 100, 100,
                                                    100, 100, 100, 100, 100,
                                                    100, 100, 100, 100, 100},
                                                    new int[] {100, 100, 100, 100, 100,
                                                    100, 100, 100, 100, 100,
                                                    100, 100, 100, 100, 100}
  };

  private static final int[][] VIRTUAL_OFFSETS = {new int[] {312, 631, 951, 1271, 1591,
                                                              1911, 2232, 2551, 2873, 3193,
                                                              3515, 3835, 4157, 4477, 4798},
                                                              new int[] {5116, 5338, 5567, 5790, 6019,
                                                              6241, 6463, 6693, 6922, 7145,
                                                              7368, 7592, 7815, 8038, 8261}
  };

  private static final int[][] VIRTUAL_OFFSET_ENDS = {new int[] {631, 951, 1271, 1591, 1911,
                                                              2232, 2551, 2873, 3193, 3515,
                                                              3835, 4157, 4477, 4798, 5116},
                                                              new int[] {5338, 5567, 5790, 6019, 6241,
                                                              6463, 6693, 6922, 7145, 7368,
                                                              7592, 7815, 8038, 8261, 100139008}
  };

  private static final int[] BINS = {5561, 4681};

  public void testSomeMethod() throws IOException {
    try (InputStream is = Resources.getResourceAsStream("com/rtg/sam/resources/mixed.sam.gz")) {
      final SamPositionReader spr = new SamPositionReader(new BlockCompressedLineReader(new BlockCompressedInputStream(is)), 0);
      try {
        int ref = 0;
        int i = 0;
        while (spr.hasNext()) {
          spr.next();
          if (i >= ENTRIES[ref]) {
            i = 0;
            ref++;
          }
          assertEquals(EXP_REF_NAME[ref], spr.getReferenceName());
          assertEquals(ref, spr.getReferenceId());
          assertEquals(START[ref][i], spr.getStartPosition());
          assertEquals(LENGTH[ref][i], spr.getLengthOnReference());
          assertEquals(BINS[ref], spr.getBinNum());
          assertEquals(VIRTUAL_OFFSETS[ref][i], spr.getVirtualOffset());
          assertEquals(VIRTUAL_OFFSET_ENDS[ref][i], spr.getNextVirtualOffset());
          assertTrue(spr.hasReference());
          assertTrue(spr.hasCoordinates());
          assertFalse(spr.isUnmapped());
          i++;
        }
      } finally {
        spr.close();
      }
    }
  }
}
