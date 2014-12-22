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

package com.rtg.util.arithcode;

import java.util.Random;

import com.rtg.util.array.byteindex.ByteChunks;

import junit.framework.TestCase;

/**
 * Test of <code>ArithEncoder</code> and <code>ArithDecoder</code>.
 */
public class ArithTest extends TestCase {

  //simple test just one block
  public void testSimple() {
    final ByteChunks bc = new ByteChunks(0);
    final OutputBytes ob = new OutputBytes(bc);
    final ArithEncoder en = new ArithEncoder(ob);

    assertEquals(0, en.endBlock());
    en.encode(0, 1, 5);
    en.encode(3, 5, 5);
    en.encode(0, 6, 7);
    final long eb1 = en.endBlock();
    assertTrue(0 < eb1);
    en.close();
    assertEquals(eb1, en.endBlock());

    final ArithDecoder de = new ArithDecoder(new InputBytes(bc, 0, eb1));

    final int c1 = de.getCurrentSymbolCount(5);
    assertTrue(0 <= c1 && c1 < 1);
    de.removeSymbolFromStream(0, 1, 5);

    final int c2 = de.getCurrentSymbolCount(5);
    assertTrue(3 <= c2 && c2 < 5);
    de.removeSymbolFromStream(3, 5, 5);

    final int c3 = de.getCurrentSymbolCount(7);
    assertTrue(0 <= c3 && c2 < 6);
    de.removeSymbolFromStream(0, 6, 7);
  }

  //multiple variable length blocks
  public void testBig() {
    randomTest(3, 1, 142, 1132, new StaticModel(new int[] {3, 2, 7}));
    randomTest(13, 10, 42, 132, new UniformModel(13));
    randomTest(3, 3, 142, 1132, new StaticModel(new int[] {3, 2, 7}));
  }

  /**
   *
   * @param range of individual symbols.
   * @param blocks number of different blocks.
   * @param seedi for lengths of blocks.
   * @param seedj for symbols being generated.
   * @param am model capable of handling range.
   */
  private void randomTest(final int range, final int blocks, final int seedi, final int seedj, final ArithCodeModel am) {
    final Random randi = new Random(seedi);
    final Random randj = new Random(seedj);
    final int[][] ra = new int[blocks][];
    for (int i = 0; i < ra.length; i++) {
      ra[i] = new int[randi.nextInt(50)];
      for (int j = 0; j < ra[i].length; j++) {
        ra[i][j] = randj.nextInt(range);
      }
    }

    final ByteChunks bc = new ByteChunks(0);
    final ArithEncoder en = new ArithEncoder(new OutputBytes(bc));
    final long[] positions = new long[ra.length + 1];
    positions[0] = en.endBlock();
    for (int i = 0; i < ra.length; i++) {
      for (int j = 0; j < ra[i].length; j++) {
        am.encode(en, ra[i][j]);
      }
      positions[i + 1] = en.endBlock();
    }
    //    for (int i = 0; i < ra.length; i++) {
    //      System.err.println("[" + i + "] " + positions[i] + "  " + ra[i].length);
    //    }
    //    System.err.println("[" + ra.length + "] " + positions[ra.length]);

    en.close();
    for (int i = 0; i < ra.length; i++) {
      final InputBytes ib = new InputBytes(bc, positions[i], positions[i + 1]);
      final ArithDecoder de = new ArithDecoder(ib);
      for (int j = 0; j < ra[i].length; j++) {
        final int sym = am.decode(de);
        assertEquals("i=" + i + " j=" + j, ra[i][j], sym);
      }
    }
  }

}
