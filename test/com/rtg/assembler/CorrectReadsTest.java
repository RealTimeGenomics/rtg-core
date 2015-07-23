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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class CorrectReadsTest extends TestCase {
  public void testCorrection() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final String[] left = new String[10];
      final String[] right = new String[10];
      for (int i = 0; i < left.length; i++) {
        left[i] = "CCCAGGAGAGG";
        right[i] = "AACGGGGGGTTTTAT";
      }
      left[0] = "CCCAGGTGAGG";
      right[0] = "AACGTGGGGTTTTAT";
      final File in = new File(tmpDir, "in");
      ReaderTestUtils.createPairedReaderDNA(ReaderTestUtils.fasta(left), ReaderTestUtils.fasta(right), in, new SdfId());

      final File out = new File(tmpDir, "out");
      CorrectReads.correct(in, out, 6, 2);
//      System.err.println(Arrays.toString(out.listFiles()));
      try (SequencesReader leftReader = SequencesReaderFactory.createDefaultSequencesReader(new File(out, "left"))) {
        assertEquals("CCCAGGAGAGG", DnaUtils.bytesToSequenceIncCG(leftReader.read(0)));
      }

      try (SequencesReader rightReader = SequencesReaderFactory.createDefaultSequencesReader(new File(out, "right"))) {
        assertEquals("AACGGGGGGTTTTAT", DnaUtils.bytesToSequenceIncCG(rightReader.read(0)));
      }


    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testCorrectionSingleEnd() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final String[] left = new String[10];
      for (int i = 0; i < left.length; i++) {
        left[i] = "CCCAGGAGAGG";
      }
      left[0] = "CCCAGGTGAGG";
      final File in = ReaderTestUtils.getDNADir(ReaderTestUtils.fasta(left), new File(tmpDir, "input"));

      final File out = new File(tmpDir, "out");
      CorrectReads.correct(in, out, 6, 2);
//      System.err.println(Arrays.toString(out.listFiles()));
      try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(new File(out, "0"))) {
        for (long seq = 0; seq < reader.numberSequences(); seq++) {
          assertEquals("CCCAGGAGAGG", DnaUtils.bytesToSequenceIncCG(reader.read(seq)));
        }
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
