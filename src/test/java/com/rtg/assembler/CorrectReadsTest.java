/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      for (int i = 0; i < left.length; ++i) {
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
      for (int i = 0; i < left.length; ++i) {
        left[i] = "CCCAGGAGAGG";
      }
      left[0] = "CCCAGGTGAGG";
      final File in = ReaderTestUtils.getDNADir(ReaderTestUtils.fasta(left), new File(tmpDir, "input"));

      final File out = new File(tmpDir, "out");
      CorrectReads.correct(in, out, 6, 2);
//      System.err.println(Arrays.toString(out.listFiles()));
      try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(new File(out, "0"))) {
        for (long seq = 0; seq < reader.numberSequences(); ++seq) {
          assertEquals("CCCAGGAGAGG", DnaUtils.bytesToSequenceIncCG(reader.read(seq)));
        }
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
