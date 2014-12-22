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

import static com.rtg.mode.SequenceMode.BIDIRECTIONAL;
import static com.rtg.mode.SequenceMode.PROTEIN;
import static com.rtg.mode.SequenceMode.TRANSLATED;
import static com.rtg.mode.SequenceMode.UNIDIRECTIONAL;

import java.io.ObjectStreamException;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class SequenceModeTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.SequenceMode()}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(SequenceMode.class, "[PROTEIN, BIDIRECTIONAL, UNIDIRECTIONAL, TRANSLATED]");
  }

  public void testReadResolve() throws ObjectStreamException {
    for (SequenceMode t : SequenceMode.values()) {
      assertEquals(t, t.readResolve());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode()}.
   */
  public final void testHashEquals() {
    TestUtils.equalsTest(SequenceMode.values());
    final SequenceMode[][] vg = new SequenceMode[SequenceMode.values().length][];
    for (int i = 0; i < vg.length; i++) {
      vg[i] = new SequenceMode[1];
    }
    for (final Object sm : SequenceMode.values()) {
      vg[((SequenceMode) sm).ordinal()][0] = (SequenceMode) sm;
    }
    TestUtils.equalsHashTest(vg);
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode()}.
   */
  public final void testHash() {
    for (final Object sm : SequenceMode.values()) {
      assertEquals(((SequenceMode) sm).ordinal() + 1, sm.hashCode());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#type()}.
   */
  public final void testType() {
    assertEquals(SequenceType.PROTEIN, PROTEIN.type());
    assertEquals(SequenceType.DNA, BIDIRECTIONAL.type());
    assertEquals(SequenceType.DNA, UNIDIRECTIONAL.type());
    assertEquals(SequenceType.DNA, TRANSLATED.type());
    assertEquals("PROTEIN", SequenceMode.PROTEIN.name());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#type()}.
   */
  public final void testCodeType() {
    assertEquals(SequenceType.PROTEIN, PROTEIN.codeType());
    assertEquals(SequenceType.DNA, BIDIRECTIONAL.codeType());
    assertEquals(SequenceType.DNA, UNIDIRECTIONAL.codeType());
    assertEquals(SequenceType.PROTEIN, TRANSLATED.codeType());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#numberFrames()}.
   */
  public final void testNumberFrames() {
    assertEquals(1, PROTEIN.numberFrames());
    assertEquals(2, BIDIRECTIONAL.numberFrames());
    assertEquals(1, UNIDIRECTIONAL.numberFrames());
    assertEquals(6, TRANSLATED.numberFrames());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#numberFrames()}.
   */
  public final void testCodeIncrement() {
    assertEquals(1, PROTEIN.codeIncrement());
    assertEquals(1, BIDIRECTIONAL.codeIncrement());
    assertEquals(1, UNIDIRECTIONAL.codeIncrement());
    assertEquals(3, TRANSLATED.codeIncrement());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#frameFromCode(int)}.
   */
  public final void testFrameFromCode() {
    assertEquals(ProteinFrame.PROTEIN, PROTEIN.frameFromCode(0));
    assertEquals(UnidirectionalFrame.FORWARD, UNIDIRECTIONAL.frameFromCode(0));
    assertEquals(BidirectionalFrame.FORWARD, BIDIRECTIONAL.frameFromCode(0));
    assertEquals(BidirectionalFrame.REVERSE, BIDIRECTIONAL.frameFromCode(1));
    assertEquals(TranslatedFrame.FORWARD1, TRANSLATED.frameFromCode(0));
    assertEquals(TranslatedFrame.FORWARD2, TRANSLATED.frameFromCode(1));
    assertEquals(TranslatedFrame.FORWARD3, TRANSLATED.frameFromCode(2));
    assertEquals(TranslatedFrame.REVERSE1, TRANSLATED.frameFromCode(3));
    assertEquals(TranslatedFrame.REVERSE2, TRANSLATED.frameFromCode(4));
    assertEquals(TranslatedFrame.REVERSE3, TRANSLATED.frameFromCode(5));

    try {
      assertEquals(TranslatedFrame.REVERSE3, TRANSLATED.frameFromCode(6));
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      // expected
    }

  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode#valueOf(String)}.
   */
  public final void testValueOf() {
    try {
      SequenceMode.valueOf("XX");
      fail("IllegalArgumentException expected");
    } catch (final IllegalArgumentException e) {
      assertEquals("XX", e.getMessage());

    }

    try {
      SequenceMode.valueOf("");
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }

  }

  /**
   * Test method for {@link com.rtg.mode.SequenceMode}}.
   * Checks the methods that deal with converting internal/external ids.
   */
  public final void testIntExt() {
    for (final Object sm : SequenceMode.values()) {
      final Frame[] frames = ((SequenceMode) sm).allFrames();
      for (final Frame fr : frames) {
        checkIntExt((SequenceMode) sm, fr, 0, 0);
        checkIntExt((SequenceMode) sm, fr, 25, 24);
        checkIntExt((SequenceMode) sm, fr, (Integer.MAX_VALUE / 6) - 1, 0);
      }
    }
    try {
      checkIntExt(SequenceMode.TRANSLATED, TranslatedFrame.FORWARD3, Integer.MAX_VALUE / 6, 0);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      assertEquals("Internal id out of range. id=357913941 frame=FORWARD3 frame ord=2 internal id=2147483648", e.getMessage());
    }
  }

  private void checkIntExt(final SequenceMode sm, final Frame fr, final long externalId, final long offset) {
    final int internal = sm.internalId(externalId, offset, fr);
    final Frame fra = sm.frame(internal);
    assertEquals(fr, fra);
    final long exta = sm.sequenceId(internal, offset);
    assertEquals(externalId, exta);
  }


  public void testRE() {
    try {
      SequenceMode.BIDIRECTIONAL.internalId(1, 5, BidirectionalFrame.FORWARD);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Internal id out of range. id=1 frame=FORWARD internal id=-8", e.getMessage());
    }
  }
}

