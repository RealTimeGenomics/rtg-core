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


import com.rtg.mode.DnaUtils;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceType;
import com.rtg.reader.CompressedMemorySequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class SharedProteinResourcesTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void test2() throws Exception {
    final CompressedMemorySequencesReader t = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray("actg".getBytes())}, new String[] {"seq1"}, new long[] {4}, 4, 4, SequenceType.PROTEIN);
    final CompressedMemorySequencesReader r = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray("ttttt".getBytes())}, new String[] {"seq2"}, new long[] {5}, 5, 5, SequenceType.DNA);
    final ProteinScoringMatrix m = new ProteinScoringMatrix();
    final SharedProteinResources res = new SharedProteinResources(m, t, r, true);
    assertEquals(t, res.templateReader());
    assertEquals(r, res.queryReader());
    assertEquals(m, res.proteinScoringMatrix());
    assertEquals(4, res.templateLength(0));
    assertEquals(5, res.queryLength(0));
    assertEquals(4, res.totalTemplateLength());
    assertEquals(1, res.query(0, 1).length);
    assertEquals(15, res.query(0, 1)[0]);
    assertEquals(15, res.query(0, 2)[0]);
    assertEquals(15, res.query(0, 3)[0]);
    assertEquals(13, res.query(0, -1)[0]);
    assertEquals(13, res.query(0, -2)[0]);
    assertEquals(13, res.query(0, -3)[0]);
    assertEquals("seq1", res.templateNames().name(0));
    assertEquals("seq2", res.readNames().name(0));

    assertEquals(1, res.template(0)[0]);
    assertEquals(2, res.template(0)[1]);
    assertEquals(4, res.template(0)[2]);
    assertEquals(3, res.template(0)[3]);
  }

  public void test3() throws Exception {
    final CompressedMemorySequencesReader t = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray("actg".getBytes())}, new String[] {"seq1"}, new long[] {4}, 4, 4, SequenceType.PROTEIN);
    final CompressedMemorySequencesReader r = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray("ttttt".getBytes())}, new String[] {"seq2"}, new long[] {5}, 5, 5, SequenceType.DNA);
    final ProteinScoringMatrix m = new ProteinScoringMatrix();
    final SharedProteinResources res = new SharedProteinResources(m, t, r, false);
    assertNull(res.readNames());
  }



}
