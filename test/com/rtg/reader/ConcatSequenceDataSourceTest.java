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

import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.Residue;
import com.rtg.mode.SequenceType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.InformationEvent;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests corresponding class
 */

public class ConcatSequenceDataSourceTest extends TestCase {
  private static final Residue A = DNA.A;
  private static final Residue C = DNA.C;
  private static final Residue G = DNA.G;
  private static final Residue T = DNA.T;
  private static final Residue N = DNA.N;

  public static Test suite() {
    return new TestSuite(ConcatSequenceDataSourceTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testEmpty() {
    Diagnostic.setLogStream();
    try {
      assertNotNull(new ConcatSequenceDataSource<>(null, null));
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("Cannot concatenate 0 sources", e.getMessage());
    }

    try {
      assertNotNull(new ConcatSequenceDataSource<>(new ArrayList<SequenceDataSource>(), null));
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("Cannot concatenate 0 sources", e.getMessage());
    }
  }

  public void testSingle() throws Exception {
    final byte[][] data = {{(byte) A.ordinal(), (byte) C.ordinal(), (byte) G.ordinal(), (byte) T.ordinal()}, {(byte) A.ordinal(), (byte) A.ordinal(), (byte) N.ordinal()}};
    final String[] labels = {"acgt", "aan"};
    final ArrayList<SequenceDataSource> list = new ArrayList<>();
    list.add(new ArraySequenceDataSource(data, null, labels, SequenceType.DNA));
    final String[] names = {"blah", "rah", "grah" };

    final int[] infoEvents = {0};

    final DiagnosticListener dl = new DiagnosticListener() {
      @Override
      public void close() { }

      @Override
      public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
        assertTrue(event instanceof InformationEvent);
        assertEquals("Processing \"blah\" (1 of 1)", event.getMessage());
        infoEvents[0]++;
      }
    };
    Diagnostic.addListener(dl);
    try {
      final ConcatSequenceDataSource<SequenceDataSource> ds = new ConcatSequenceDataSource<>(list, Arrays.asList(names));
      assertEquals(1, infoEvents[0]);
      assertEquals(0, ds.getWarningCount());
      assertTrue(ds.nextSequence());
      assertEquals(SequenceType.DNA, ds.type());

      final byte[] b = ds.sequenceData();
      assertEquals(A.ordinal(), b[0]);

      assertFalse(ds.nextSequence());

      new ConcatSequenceDataSource<>(list, null);
      assertEquals(1, infoEvents[0]);
    } finally {
      Diagnostic.removeListener(dl);
    }

  }

  private class ArraySequenceDataSource implements SequenceDataSource {
    final byte[][] mData;
    final byte[][] mQuality;
    final String[] mLabels;
    final SequenceType mType;
    int mSequenceIndex;
    long mMinLength = Long.MAX_VALUE;
    long mMaxLength = Long.MIN_VALUE;

    public ArraySequenceDataSource(byte[][] data, byte[][] quality, String[] labels, SequenceType type) {
      mData = Arrays.copyOf(data, data.length);
      mQuality = quality == null ? null : Arrays.copyOf(quality, quality.length);
      mLabels = Arrays.copyOf(labels, labels.length);
      mType = type;

      mSequenceIndex = 0;
    }

    @Override
    public SequenceType type() {
      return mType;
    }

    @Override
    public boolean nextSequence() {
      mSequenceIndex++;
      if (mSequenceIndex < mData.length) {
        mMinLength = Math.min(mMinLength, currentLength());
        mMaxLength = Math.max(mMaxLength, currentLength());
      }
      return mSequenceIndex < mData.length;
    }

    @Override
    public String name() {
      return mLabels[mSequenceIndex];
    }

    @Override
    public byte[] sequenceData() {
      return mData[mSequenceIndex];
    }

    @Override
    public byte[] qualityData() {
      return mQuality == null ? null : mQuality[mSequenceIndex];
    }

    @Override
    public boolean hasQualityData() {
      return mQuality != null;
    }

    @Override
    public void close() {
      // nothing to do
    }

    @Override
    public void setDusting(final boolean val) {
      // ignored
    }

    @Override
    public long getWarningCount() {
      return 0;
    }

    @Override
    public int currentLength() {
      return mData[mSequenceIndex].length;
    }

    @Override
    public long getDusted() {
      return 0;
    }

    @Override
    public long getMaxLength() {
      return mMaxLength;
    }

    @Override
    public long getMinLength() {
      return mMinLength;
    }
  }
}
