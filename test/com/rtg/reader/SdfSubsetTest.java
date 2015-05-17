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


import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.TestDirectory;

/**
 * Test for SdfSubset
 *
 */
public class SdfSubsetTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SdfSubset();
  }

  public void testHelp() {
    checkHelp("Extracts a subset of sequences from one SDF and outputs them to another SDF.",
              "output SDF",
              "SDF containing sequences",
              "id of sequence",
              "file containing sequence ids, or sequence names if --names flag is set, one per line");
  }

  public void testValidator() throws Exception {
    final String err = checkHandleFlagsErr("-i", "genome", "-o", "out");
    TestUtils.containsAll(err, "Sequences to extract must be specified");
  }

  private void checkPairedSplittiness(String... extraParams) throws Exception {
    try (final TestDirectory tempDir = new TestDirectory("sdfsubsettest")) {
      final String left = ">seq1" + StringUtils.LS + "acgt" + StringUtils.LS + ">seq2" + StringUtils.LS + "cgta";
      final String right = ">seq1" + StringUtils.LS + "tgca" + StringUtils.LS + ">seq2" + StringUtils.LS + "atgc";


      final File readLeftDir = new File(tempDir, "left");
      ReaderTestUtils.getReaderDNA(left, readLeftDir, null).close();
      final File readRightDir = new File(tempDir, "right");
      ReaderTestUtils.getReaderDNA(right, readRightDir, null).close();
      final File outDir = new File(tempDir, "outDir");


      checkMainInitOk(Utils.append(new String[] {"-i", tempDir.getPath(), "-o", outDir.getPath()}, extraParams));

      SdfId leftGuid;
      SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "left"));
      try {
        assertEquals(4, dsr.totalLength());
        assertEquals(4, dsr.maxLength());
        assertEquals(4, dsr.minLength());
        assertEquals(1, dsr.numberSequences());
        assertNotNull(dsr.getSdfId());
        leftGuid = dsr.getSdfId();
        assertTrue(dsr.hasNames());
        assertEquals("seq2", dsr.name(0));
        assertFalse(dsr.hasQualityData());
      } finally {
        dsr.close();
      }
      dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "right"));
      try {
        assertEquals(4, dsr.totalLength());
        assertEquals(4, dsr.maxLength());
        assertEquals(4, dsr.minLength());
        assertEquals(1, dsr.numberSequences());
        assertNotNull(dsr.getSdfId());
        assertEquals(leftGuid, dsr.getSdfId());
        assertTrue(dsr.hasNames());
        assertEquals("seq2", dsr.name(0));
        assertFalse(dsr.hasQualityData());
      } finally {
        dsr.close();
      }
    }

  }
  public void testPairedSplittiness() throws Exception {
    checkPairedSplittiness("1");
  }
  public void testPairedSplittinessNames() throws Exception {
    checkPairedSplittiness("-n", "seq2");
  }
  public void testPairedSplittinessRange() throws Exception {
    checkPairedSplittiness("--start-id", "1", "--end-id", "2");
  }
}
