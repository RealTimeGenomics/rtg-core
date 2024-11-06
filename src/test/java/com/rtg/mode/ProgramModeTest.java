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
package com.rtg.mode;

import static com.rtg.mode.ProgramMode.DUSTER;
import static com.rtg.mode.ProgramMode.PHYLOGENY;
import static com.rtg.mode.ProgramMode.SLIMN;
import static com.rtg.mode.ProgramMode.SLIMN_FLIP;
import static com.rtg.mode.ProgramMode.SLIMP;
import static com.rtg.mode.ProgramMode.SLIMX;
import static com.rtg.mode.ProgramMode.TSLIMN;
import static com.rtg.mode.ProgramMode.TSLIMX;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class ProgramModeTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.ProgramMode()}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(ProgramMode.class, "[SLIMN, SLIMP, SLIMX, TSLIMN, TSLIMX]");
  }

  /**
   * Test method for {@link com.rtg.mode.ProgramMode}.
   */
  public final void testMode() {
    assertEquals(SequenceMode.PROTEIN, SLIMP.queryMode());
    assertEquals(SequenceMode.PROTEIN, SLIMP.subjectMode());
    assertEquals(SequenceType.PROTEIN, SLIMP.translatedType());
    assertEquals(SLIMP, SLIMP.flip());
    assertEquals("SLIMP", SLIMP.name());

    assertEquals(SequenceMode.TRANSLATED, SLIMX.queryMode());
    assertEquals(SequenceMode.PROTEIN, SLIMX.subjectMode());
    assertEquals(SequenceType.PROTEIN, SLIMX.translatedType());
    assertEquals(TSLIMN, SLIMX.flip());

    assertEquals(SequenceMode.PROTEIN, TSLIMN.queryMode());
    assertEquals(SequenceMode.TRANSLATED, TSLIMN.subjectMode());
    assertEquals(SequenceType.PROTEIN, TSLIMN.translatedType());
    assertEquals(SLIMX, TSLIMN.flip());

    assertEquals(SequenceMode.TRANSLATED, TSLIMX.queryMode());
    assertEquals(SequenceMode.TRANSLATED, TSLIMX.subjectMode());
    assertEquals(SequenceType.PROTEIN, TSLIMX.translatedType());
    assertEquals(TSLIMX, TSLIMX.flip());

    assertEquals(SequenceMode.BIDIRECTIONAL, SLIMN.queryMode());
    assertEquals(SequenceMode.UNIDIRECTIONAL, SLIMN.subjectMode());
    assertEquals(SequenceType.DNA, SLIMN.translatedType());
    assertEquals(SLIMN_FLIP, SLIMN.flip());

    assertEquals(SequenceMode.UNIDIRECTIONAL, SLIMN_FLIP.queryMode());
    assertEquals(SequenceMode.BIDIRECTIONAL, SLIMN_FLIP.subjectMode());
    assertEquals(SequenceType.DNA, SLIMN_FLIP.translatedType());
    assertEquals(SLIMN, SLIMN_FLIP.flip());

    assertEquals(DUSTER, DUSTER.flip());
    assertEquals(5, DUSTER.ordinal());

    assertEquals(6, SLIMN_FLIP.ordinal());

    assertEquals(SequenceMode.UNIDIRECTIONAL, PHYLOGENY.queryMode());
    assertEquals(SequenceMode.UNIDIRECTIONAL, PHYLOGENY.subjectMode());
    assertEquals(SequenceType.DNA, PHYLOGENY.translatedType());
    assertEquals(PHYLOGENY, PHYLOGENY.flip());
    assertEquals(7, PHYLOGENY.ordinal());
  }

  public void testReadResolve() {
    assertEquals(DUSTER, DUSTER.readResolve());
    assertEquals(SLIMN_FLIP, SLIMN_FLIP.readResolve());
    for (Object t : ProgramMode.values()) {
      assertEquals(t, ((ProgramMode) t).readResolve());
    }
  }

}
