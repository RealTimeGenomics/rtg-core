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
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test for SdfSubset
 *
 */
public class SdfSubsetTest extends AbstractCliTest {

  public static Test suite() {
    return new TestSuite(SdfSubsetTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.runAndWait(suite());
  }

  @Override
  protected AbstractCli getCli() {
    return new SdfSubset();
  }

  public void testHelp() {
    checkHelp("Extracts a subset of sequences from one SDF and outputs them to another SDF.",
              "output SDF",
              "input SDF",
              "id of sequence",
              "file containing sequence ids, or sequence names if names flag is set, one per line");
  }

  public void testValidator() throws Exception {
    final File tempDir = FileUtils.createTempDir("readsimclitest", "checkcli");
    try {
      final File genomeDir = FileHelper.createTempDirectory();
      try {
        final File out = FileUtils.createTempDir("out", "validate");
        assertTrue(FileHelper.deleteAll(out));
        try {
          ReaderTestUtils.getReaderDNA(">seq1" + StringUtils.LS + "acgt", genomeDir, null).close();
          TestUtils.containsAll(checkHandleFlagsErr("-i", genomeDir.getAbsolutePath(), "-o", out.getPath()), "Sequence ids must be specified, either explicitly, or using --id-file");
        } finally {
          assertTrue(FileHelper.deleteAll(out));
        }
      } finally {
        assertTrue(FileHelper.deleteAll(genomeDir));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void checkPairedSplittiness(String... extraParams) throws Exception {
    final File tempDir = FileUtils.createTempDir("readsimclitest", "checkcli");
    try {
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
        assertTrue(dsr.nextSequence());
        assertNotNull(dsr.getSdfId());
        leftGuid = dsr.getSdfId();
        assertTrue(dsr.hasNames());
        assertEquals("seq2", dsr.currentName());
        assertFalse(dsr.hasQualityData());
        assertFalse(dsr.nextSequence());
      } finally {
        dsr.close();
      }
      dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "right"));
      try {
        assertEquals(4, dsr.totalLength());
        assertEquals(4, dsr.maxLength());
        assertEquals(4, dsr.minLength());
        assertTrue(dsr.nextSequence());
        assertNotNull(dsr.getSdfId());
        assertEquals(leftGuid, dsr.getSdfId());
        assertTrue(dsr.hasNames());
        assertEquals("seq2", dsr.currentName());
        assertFalse(dsr.hasQualityData());
        assertFalse(dsr.nextSequence());
      } finally {
        dsr.close();
      }

    } finally {
      FileHelper.deleteAll(tempDir);
    }

  }
  public void testPairedSplittiness() throws Exception {
    checkPairedSplittiness("1");
  }
  public void testPairedSplittinessNames() throws Exception {
    checkPairedSplittiness("-n", "seq2");
  }
}
