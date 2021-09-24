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
