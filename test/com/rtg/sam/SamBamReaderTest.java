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
package com.rtg.sam;

import java.util.ArrayList;

import com.rtg.sam.SamBamReader.SequenceInfo;

import junit.framework.TestCase;


/**
 * Test AbstractSamBamReader.
 *
 */
public class SamBamReaderTest extends TestCase {

  private ConcreteSamBamReader mReader;

  @Override
  public void setUp() {
    mReader = new ConcreteSamBamReader();
  }

  @Override
  public void tearDown() {
    mReader = null;
  }

  class ConcreteSamBamReader extends SamBamReader {

    ConcreteSamBamReader() {
      mSeqInfo = new ArrayList<>();
      mSeqInfo.add(new SequenceInfo("seq1", 42));
      mSeqInfo.add(new SequenceInfo("seq2", 43));
    }

    @Override
    public Object getAttributeValue(String tag) {
      return null;
    }

    @Override
    public String[] getAttributeTags() {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getFieldNumFromTag(String tag) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getField(int fieldNum) {
      return null;
    }

    @Override
    public byte[] getFieldBytes(int fieldNum) {
      return null;
    }

    @Override
    public String getHeaderLines() {
      return "@header1";
    }

    @Override
    public String getFullHeader() {
      return "@header1";
    }

    @Override
    public int getIntAttribute(String tag) {
      return 0;
    }

    @Override
    public int getIntField(int fieldNum) {
      return 0;
    }

    @Override
    public int getNumFields() {
      return 0;
    }

    @Override
    public boolean hasAttribute(String tag) {
      return false;
    }

    @Override
    public boolean hasNext() {
      return false;
    }

    @Override
    public boolean isSam() {
      return false;
    }

    @Override
    public SamBamRecord next() {
      return null;
    }

    @Override
    public void remove() {
    }

    @Override
    public void close() {
    }

    @Override
    public char getAttributeType(String tag) {
      return (char) 0;
    }

  }

  public void testGetHeaderLines() {
    assertEquals("@header1", mReader.getHeaderLines());
  }

  public void testIsSam() {
    assertFalse(mReader.isSam());
  }

  public void testNumReferences() {
    assertEquals(2, mReader.numReferences());
  }

  public void testLookupSequence() {
    assertEquals("seq1", mReader.lookupSequence(0).getName());
    assertEquals("seq2", mReader.lookupSequence(1).getName());
    try {
      mReader.lookupSequence(2);
      fail("IndexOutOfBoundsException expected");
    } catch (final IndexOutOfBoundsException e) {
      // good.
    }
  }

  public void testSequenceDictionary() {
    assertEquals(2, mReader.sequenceDictionary().size());
  }

  public void testIterator() {
    assertEquals(mReader, mReader.iterator());
  }

  public void testSequenceInfo() {
    final SequenceInfo blah = new SequenceInfo("blah", 123);
    assertEquals("blah", blah.getName());
    assertEquals(123, blah.getLength());
  }
}
