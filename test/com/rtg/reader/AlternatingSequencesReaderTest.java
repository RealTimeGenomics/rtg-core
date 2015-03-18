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

import java.io.IOException;

import com.rtg.mode.SequenceType;

import junit.framework.TestCase;

/**
 */
public class AlternatingSequencesReaderTest extends TestCase {

  public void testExceptions1() throws IllegalStateException, IOException {
    final SequencesReader first = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final SequencesReader second = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final SequencesIterator reader = new AlternatingSequencesReader(first, second).iterator();
    try {
      reader.currentLength();
      fail();
    } catch (IllegalStateException e) {
      //expected
    }
    try {
      reader.currentSequenceId();
      fail();
    } catch (IllegalStateException e) {
      //expected
    }
  }

  public void testExceptions3() throws IllegalStateException, IOException {
    final SequencesReader first = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final SequencesReader second = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final AlternatingSequencesReader reader = new AlternatingSequencesReader(first, second);
    try {
      reader.dataChecksum();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.qualityChecksum();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.nameChecksum();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.copy();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.path();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.getArm();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.getPrereadType();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.getSdfId();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.hasHistogram();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.hasQualityData();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.histogram();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.lengthBetween(0, 0);
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.longestNBlock();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.maxLength();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.minLength();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported yet.", e.getMessage());
    }
    reader.close();
  }

  public void testExceptions2() throws IllegalStateException, IOException {
    final SequencesReader first = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final SequencesReader second = new MockSequencesReader(SequenceType.DNA, 2L, 10L);
    final AlternatingSequencesReader reader = new AlternatingSequencesReader(first, second);
    try {
      reader.nBlockCount();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.names();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.posHistogram();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.positionQualityAverage();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.residueCounts();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.sdfVersion();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.sequenceLengths(0, 0);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.type();
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    reader.close();
  }

  private static class FakeReader extends MockSequencesReader {
    private final int mReadLength;

    public FakeReader(SequenceType sequenceType, long numberSquences, long length, int readLength) {
      super(sequenceType, numberSquences, length);
      mReadLength = readLength;
    }
    @Override
    public int length(long index) {
      return mReadLength + (int) index;
    }
  }

  public void test() throws IllegalStateException, IOException {
    final SequencesReader first = new FakeReader(SequenceType.DNA, 2L, 13L, 6);
    final SequencesReader second = new FakeReader(SequenceType.DNA, 2L, 11L, 5);
    final AlternatingSequencesReader reader = new AlternatingSequencesReader(first, second);
    assertEquals(4, reader.numberSequences());
    assertEquals(24, reader.totalLength());
    final SequencesIterator it = reader.iterator();
    assertTrue(it.nextSequence());
    assertEquals(0, it.currentSequenceId());
    assertEquals(6, it.currentLength());
    assertTrue(it.nextSequence());
    assertEquals(1, it.currentSequenceId());
    assertEquals(5, it.currentLength());
    assertTrue(it.nextSequence());
    assertEquals(2, it.currentSequenceId());
    assertEquals(7, it.currentLength());
    assertTrue(it.nextSequence());
    assertEquals(3, it.currentSequenceId());
    assertEquals(6, it.currentLength());
    assertFalse(it.nextSequence());
    reader.close();
  }
}
