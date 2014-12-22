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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceType;
import com.rtg.reader.FastqSequenceDataSource.FastQScoreType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ReadHelperTest extends TestCase {

  /**
   */
  public ReadHelperTest(String name) {
    super(name);
  }

  public void testHandling() throws Exception {
    final String seqString = "acgtcacgtcacgtcacgtcacgtcacgtcacgtc";
    Diagnostic.setLogStream();
    assertNull(ReadHelper.getRead(null, 0));
    assertNull(ReadHelper.getQual(null, 0));
    CompressedMemorySequencesReader msr = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray(seqString.getBytes())}, new String[] {"seq1"}, new long[] {35}, 35, 35, SequenceType.DNA) {
      private boolean mIsLeft = true;

      @Override
      public PrereadArm getArm() {
        try {
          return mIsLeft ? PrereadArm.LEFT : PrereadArm.RIGHT;
        } finally {
          mIsLeft = false;
        }
      }

      @Override
      public PrereadType getPrereadType() {
        return PrereadType.UNKNOWN;
      }

      @Override
      public boolean equals(Object o) {
        return o == this;
      }

      @Override
      public int hashCode() {
        return 1;
      }



    };
    try {
      assertTrue(Arrays.equals(DnaUtils.encodeArray("acgtcacgtcacgtcacgtcacgtcacgtcacgtc".getBytes()), ReadHelper.getRead(msr, 0)));
      assertTrue(Arrays.equals(DnaUtils.encodeArray("acgtcacgtcacgtcacgtcacgtcacgtcacgtc".getBytes()), ReadHelper.getRead(msr, 0)));
      assertNull(ReadHelper.getQual(msr, 0));
    } finally {
      msr.close();
    }
    msr = new CompressedMemorySequencesReader(new byte[][] {DnaUtils.encodeArray(seqString.getBytes())}, new String[] {"seq1"}, new long[] {35}, 35, 35, SequenceType.DNA) {
      private boolean mIsLeft = true;

      @Override
      public PrereadArm getArm() {
        try {
          return mIsLeft ? PrereadArm.LEFT : PrereadArm.RIGHT;
        } finally {
          mIsLeft = false;
        }
      }

      @Override
      public PrereadType getPrereadType() {
        return PrereadType.UNKNOWN;
      }

      @Override
      public boolean equals(Object o) {
        return o == this;
      }

      @Override
      public int hashCode() {
        return 1;
      }
    };
    try {
      assertTrue(Arrays.equals(DnaUtils.encodeArray("acgtcacgtcacgtcacgtcacgtcacgtcacgtc".getBytes()), ReadHelper.getRead(msr, 0)));
    } finally {
      msr.close();
    }
  }

  public void testCGQuality() throws Exception {
    Diagnostic.setLogStream();
    final File temp = FileUtils.createTempDir("cgblah", "qual");
    try {
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(new ByteArrayInputStream(("@testQuality\n"
          + "actgcatc\n"
          + "+\n"
          + "!<><##!<").getBytes()));
      final FastqSequenceDataSource fq = new FastqSequenceDataSource(al, FastQScoreType.PHRED);
      final SequencesWriter sw = new SequencesWriter(fq, temp, 20, PrereadType.CG, false);
      sw.processSequences();
      try (SequencesReader msr = CompressedMemorySequencesReader.createSequencesReader(temp, true, false, LongRange.NONE)) {
        assertTrue(Arrays.equals(DnaUtils.encodeArray("actgcatc".getBytes()), ReadHelper.getRead(msr, 0)));
        assertTrue(Arrays.equals(new byte[]{0, 27, 29, 27, 2, 2, 0, 27}, ReadHelper.getQual(msr, 0)));
      }
    } finally {
      FileHelper.deleteAll(temp);
    }
  }
}
