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

package com.rtg.simulation.reads;

import java.io.IOException;

import com.rtg.reader.SdfId;
import com.rtg.util.PortableRandom;

import junit.framework.TestCase;

/**
 */
public class UnknownBaseReadWriterTest extends TestCase {

  private static class MockReadWriter implements ReadWriter {

    int mFrom = 0;
    byte[] mLastData;
    @Override
    public void close() {
    }

    @Override
    public void identifyTemplateSet(SdfId... templateIds) {
    }

    @Override
    public void identifyOriginalReference(SdfId referenceId) {
    }

    @Override
    public void writeRead(String name, byte[] data, byte[] qual, int length) {
      mLastData = new byte[length];
      System.arraycopy(data, 0, mLastData, 0, length);
      mFrom = 0;
    }

    @Override
    public void writeLeftRead(String name, byte[] data, byte[] qual, int length) {
      writeRead(name, data, qual, length);
      mFrom = 1;
    }

    @Override
    public void writeRightRead(String name, byte[] data, byte[] qual, int length) {
      writeRead(name, data, qual, length);
      mFrom = 2;
    }

    @Override
    public int readsWritten() {
      return 79;
    }
  }
  static class TestRandom  extends PortableRandom {
      int mI = 0;
      double[] mResults = {0.7, 0.1, 0.8, 0.9, 0.5, 0.2, 0.6, 0.1, 0.3, 0.4, 0.0};
      @Override
      public double nextDouble() {
        return mResults[mI++ % mResults.length];
      }
  }

  public void test() throws IOException {
    MockReadWriter internal = new MockReadWriter();
    ReadWriter rw = new UnknownBaseReadWriter(internal, 0.3, new TestRandom());
    byte[] in = {1, 2, 3, 4, 5, 6, 7, 8};
    rw.writeLeftRead("asdf", in, new byte[8], 8);
    assertEquals(1, internal.mFrom);
    byte[] expected = {1, 0, 3, 4, 5, 0, 7, 0};
    for (int i = 0; i < internal.mLastData.length; i++) {
      assertEquals(expected[i], internal.mLastData[i]);
    }
    rw.writeRightRead("asdf", in, new byte[8], 8);
    assertEquals(2, internal.mFrom);
    expected = new byte[] {1, 2, 0, 4, 0, 6, 7, 8};
    for (int i = 0; i < internal.mLastData.length; i++) {
      assertEquals(expected[i], internal.mLastData[i]);
    }
    in = new byte[] {1, 2, 3, 4, 5, 6, 7, 8, 9};
    rw.writeRightRead("asdf", in, new byte[8], 9);
    assertEquals(2, internal.mFrom);
    expected = new byte[] {0, 2, 0, 4, 5, 0, 7, 0, 9};
    for (int i = 0; i < internal.mLastData.length; i++) {
      assertEquals(expected[i], internal.mLastData[i]);
    }
  }

  public void test2() throws IOException {
    MockReadWriter internal = new MockReadWriter();
    ReadWriter rw = new UnknownBaseReadWriter(internal, 0.2, new TestRandom());
    byte[] in = {1, 2, 3, 4, 5, 6, 7, 8};
    rw.writeLeftRead("asdf", in, new byte[8], 8);
    assertEquals(1, internal.mFrom);
    byte[] expected = {1, 0, 3, 4, 5, 6, 7, 0};
    for (int i = 0; i < internal.mLastData.length; i++) {
      assertEquals(expected[i], internal.mLastData[i]);
    }
    assertEquals(79, rw.readsWritten());
  }

}
