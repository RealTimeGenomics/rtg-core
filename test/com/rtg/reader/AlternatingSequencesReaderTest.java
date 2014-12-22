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
    final AlternatingSequencesReader reader = new AlternatingSequencesReader(first, second);
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
    reader.close();
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
      reader.currentName();
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
      reader.length(0);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.read(0, null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.read(0, null, 0, 0);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.readCurrent(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.readCurrent(null, 0, 0);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.readCurrentQuality(null);
      fail();
    } catch (UnsupportedOperationException e) {
      //expected
      assertEquals("Not supported yet.", e.getMessage());
    }
    try {
      reader.readQuality(0, null);
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
    public int currentLength() {
      return mReadLength + (int) currentSequenceId();
    }
  }

  public void test() throws IllegalStateException, IOException {
    final SequencesReader first = new FakeReader(SequenceType.DNA, 2L, 13L, 6);
    final SequencesReader second = new FakeReader(SequenceType.DNA, 2L, 11L, 5);
    final AlternatingSequencesReader reader = new AlternatingSequencesReader(first, second);
    assertEquals(4, reader.numberSequences());
    assertEquals(24, reader.totalLength());
    assertTrue(reader.nextSequence());
    assertEquals(0, reader.currentSequenceId());
    assertEquals(6, reader.currentLength());
    assertTrue(reader.nextSequence());
    assertEquals(1, reader.currentSequenceId());
    assertEquals(5, reader.currentLength());
    assertTrue(reader.nextSequence());
    assertEquals(2, reader.currentSequenceId());
    assertEquals(7, reader.currentLength());
    assertTrue(reader.nextSequence());
    assertEquals(3, reader.currentSequenceId());
    assertEquals(6, reader.currentLength());
    assertFalse(reader.nextSequence());
    reader.close();
  }
}
