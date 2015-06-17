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
package com.rtg.mode;

import static com.rtg.mode.ProgramMode.DUSTER;
import static com.rtg.mode.ProgramMode.PHYLOGENY;
import static com.rtg.mode.ProgramMode.SLIMN;
import static com.rtg.mode.ProgramMode.SLIMN_FLIP;
import static com.rtg.mode.ProgramMode.SLIMP;
import static com.rtg.mode.ProgramMode.SLIMX;
import static com.rtg.mode.ProgramMode.TSLIMN;
import static com.rtg.mode.ProgramMode.TSLIMX;

import java.io.ObjectStreamException;

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

  public void testReadResolve() throws ObjectStreamException {
    assertEquals(DUSTER, DUSTER.readResolve());
    assertEquals(SLIMN_FLIP, SLIMN_FLIP.readResolve());
    for (Object t : ProgramMode.values()) {
      assertEquals(t, ((ProgramMode) t).readResolve());
    }
  }

}
