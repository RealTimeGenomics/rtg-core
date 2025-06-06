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
package com.rtg.protein;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.BidirectionalEditDistance;
import com.rtg.alignment.EditDistanceFactory;
import com.rtg.mode.Frame;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.TranslatedFrame;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ProteinAlignmentResultTest extends TestCase {

  /**
   */
  public ProteinAlignmentResultTest(String name) {
    super(name);
  }

  public void testCreation() throws IOException, InvalidParamsException {
    final File temp = FileUtils.createTempDir("prot", "alignmetresult");
    try {
      Diagnostic.setLogStream();
      final File template = new File(temp, "template");
      ReaderTestUtils.getReaderProtein(ProteinOutputProcessorTest.TEMPLATE_FASTA, template);
      final SequencesReader tem = SequencesReaderFactory.createMemorySequencesReader(template, true, LongRange.NONE);
      final File readsFile = new File(temp, "reads");
      ReaderTestUtils.getReaderDNA(ProteinOutputProcessorTest.READS_FASTA_PERFECT, readsFile, null);
      final SequencesReader reads = SequencesReaderFactory.createMemorySequencesReader(readsFile, true, LongRange.NONE);
      final BidirectionalEditDistance ed = EditDistanceFactory.createProteinEditDistance(new ProteinScoringMatrix());

      final SharedProteinResources resx = new SharedProteinResources(new ProteinScoringMatrix(), tem, reads, false);
      final byte[] t = tem.read(0);

      final byte[] b = reads.read(0);
      final int rlen = b.length;

      final int genomeFrame = 0;
      final Frame frames = TranslatedFrame.FORWARD1;
      final int plen = (rlen - Math.abs(genomeFrame) + 1) / 3;
      final byte[] p = new byte[plen];
      for (int j = 0, i = 0; j < plen; ++j, i += 3) {
        p[j] = frames.code(b, rlen, i);
      }
      final int[] r = ed.calculateEditDistance(p, plen, t, 0, false, Integer.MAX_VALUE, MaxShiftUtils.calculateDefaultMaxShift(plen), true);
      final ProteinAlignmentResult res = new ProteinAlignmentResult(resx, 0, 0, r, 0, true);
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        res.write(bos);
      } finally {
        bos.close();
      }
      assertEquals("templateName\t+1\t0\t1\t12\t22\t1\t36\t36\tkwrknrkskknq\tkwrknrkskknq\tkwrknrkskknq\t12\t100\t12\t100\t0\t-67\t30.4\t1.5e-8", bos.toString());

    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  public void testCompareTo() {

    final int[] actions = new int[ActionsHelper.ACTIONS_START_INDEX + 1];
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 100;
    final ProteinAlignmentResult r1 = new ProteinAlignmentResult(null, 10, 10, actions, 0, true);

    // Same result
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 10, 10, actions, 0, true)) == 0);

    // Ordering of template id
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 11, 10, actions, 0, true)) < 0);
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 9, 10, actions, 0, true)) > 0);

    // Ordering of start position
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 101;
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 10, 10, actions, 0, true)) < 0);
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 99;
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 10, 10, actions, 0, true)) > 0);

    // Ordering of read id
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 100;
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 10, 11, actions, 0, true)) < 0);
    assertTrue(r1.compareTo(new ProteinAlignmentResult(null, 10, 9, actions, 0, true)) > 0);
  }

  public void testEqual() {
    final int[] actions = new int[ActionsHelper.ACTIONS_START_INDEX + 1];
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 100;
    final ProteinAlignmentResult r1 = new ProteinAlignmentResult(null, 10, 10, actions, 0, true);

    assertTrue(r1.equals(new ProteinAlignmentResult(null, 10, 10, actions, 0, true)));
    assertTrue(r1.equals(r1));

    assertFalse(r1.equals(null));
    assertFalse(r1.equals(new ProteinAlignmentResult(null, 10, 11, actions, 0, true)));
    final int[] actions2 = new int[ActionsHelper.ACTIONS_START_INDEX + 1];
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 101;
    assertFalse(r1.equals(new ProteinAlignmentResult(null, 10, 10, actions2, 0, true)));
  }

  public void testHashCode() {
    final int[] actions = new int[ActionsHelper.ACTIONS_START_INDEX + 1];
    actions[ActionsHelper.TEMPLATE_START_INDEX] = 100;
    final ProteinAlignmentResult r1 = new ProteinAlignmentResult(null, 10, 10, actions, 0, true);
    assertEquals(new ProteinAlignmentResult(null, 10, 10, actions, 0, true).hashCode(), r1.hashCode());
    assertFalse(new ProteinAlignmentResult(null, 10, 11, actions, 0, true).hashCode() == r1.hashCode());
  }

  public void testReadPositions() {
    assertEquals(1, ProteinAlignmentResult.readNtStart(1, 6));
    assertEquals(2, ProteinAlignmentResult.readNtStart(2, 6));
    assertEquals(3, ProteinAlignmentResult.readNtStart(3, 6));
    assertEquals(1, ProteinAlignmentResult.readNtStart(-1, 6));
    assertEquals(3, ProteinAlignmentResult.readNtStart(-2, 6));
    assertEquals(2, ProteinAlignmentResult.readNtStart(-3, 6));
    assertEquals(6, ProteinAlignmentResult.readNtEnd(1, 6));
    assertEquals(4, ProteinAlignmentResult.readNtEnd(2, 6));
    assertEquals(5, ProteinAlignmentResult.readNtEnd(3, 6));

    assertEquals(1, ProteinAlignmentResult.readNtStart(1, 7));
    assertEquals(2, ProteinAlignmentResult.readNtStart(2, 7));
    assertEquals(3, ProteinAlignmentResult.readNtStart(3, 7));
    assertEquals(2, ProteinAlignmentResult.readNtStart(-1, 7));
    assertEquals(1, ProteinAlignmentResult.readNtStart(-2, 7));
    assertEquals(3, ProteinAlignmentResult.readNtStart(-3, 7));
    assertEquals(6, ProteinAlignmentResult.readNtEnd(1, 7));
    assertEquals(7, ProteinAlignmentResult.readNtEnd(2, 7));
    assertEquals(5, ProteinAlignmentResult.readNtEnd(3, 7));

    assertEquals(1, ProteinAlignmentResult.readNtStart(1, 8));
    assertEquals(2, ProteinAlignmentResult.readNtStart(2, 8));
    assertEquals(3, ProteinAlignmentResult.readNtStart(3, 8));
    assertEquals(3, ProteinAlignmentResult.readNtStart(-1, 8));
    assertEquals(2, ProteinAlignmentResult.readNtStart(-2, 8));
    assertEquals(1, ProteinAlignmentResult.readNtStart(-3, 8));
    assertEquals(6, ProteinAlignmentResult.readNtEnd(1, 8));
    assertEquals(7, ProteinAlignmentResult.readNtEnd(2, 8));
    assertEquals(8, ProteinAlignmentResult.readNtEnd(3, 8));

    assertEquals(1, ProteinAlignmentResult.readNtStart(1, 4));
    assertEquals(2, ProteinAlignmentResult.readNtStart(2, 4));
    assertEquals(3, ProteinAlignmentResult.readNtStart(3, 4));
    assertEquals(2, ProteinAlignmentResult.readNtStart(-1, 4));
    assertEquals(1, ProteinAlignmentResult.readNtStart(-2, 4));
    assertEquals(0, ProteinAlignmentResult.readNtStart(-3, 4)); // off end, can't actually happen
    assertEquals(3, ProteinAlignmentResult.readNtEnd(1, 4));
    assertEquals(4, ProteinAlignmentResult.readNtEnd(2, 4));
    assertEquals(2, ProteinAlignmentResult.readNtEnd(3, 4)); // correct given above
  }

  public void testWritingBitScore() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ProteinAlignmentResult.writeBitScore(bos, 0.0);
    assertEquals("0.0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, -1);
    assertEquals("-1.0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, 1e-180);
    assertEquals("0.0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, 1.2);
    assertEquals("1.2", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, 1234);
    assertEquals("1234.0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, -123456.7890);
    assertEquals("-123456.8", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeBitScore(bos, 123456.7890);
    assertEquals("123456.8", bos.toString());

  }

  public void testWritingEScore() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ProteinAlignmentResult.writeEScore(bos, 0.0);
    assertEquals("0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeEScore(bos, 1e-180);
    assertEquals("1.0e-180", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeEScore(bos, 1.2);
    assertEquals("1.2e0", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeEScore(bos, 1234);
    assertEquals("1.2e3", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeEScore(bos, 123456.7890);
    assertEquals("1.2e5", bos.toString());
    bos.reset();
    ProteinAlignmentResult.writeEScore(bos, -123456.7890);
    assertEquals("0", bos.toString());

  }


}
