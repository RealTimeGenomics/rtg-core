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
package com.rtg.variant.cnv;

import java.io.File;
import java.io.FileOutputStream;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class NBlockDetectorTest extends TestCase {

  private static final String SEQUENCES = ">seq1" + StringUtils.LS
          + "aaaannaaca"
          + "nnnnnnnnnn"
          + "nacacggttc"
          + "ancagttacn"
          + "nncagtcagc"
          + "atnn" + StringUtils.LS

          + ">seq2" + StringUtils.LS
          + "nnn" + StringUtils.LS

          + ">seq3" + StringUtils.LS
          + "nnnactacga"
          + "gcatgact" + StringUtils.LS

          + ">seq4" + StringUtils.LS
          + "acgatcagtc"
          + "agctnnn" + StringUtils.LS;

  private static final String EXPECTED = ""
          + "seq1\t11\t21" + StringUtils.LS
          + "seq1\t40\t42" + StringUtils.LS
          + "seq2\t1\t3" + StringUtils.LS
          + "seq3\t1\t3" + StringUtils.LS
          + "seq4\t15\t17" + StringUtils.LS;

  public void testDetectNs() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("nblockdetect", "test");
    try {
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getReaderDNA(SEQUENCES, sdf, null).close();
      final File output = new File(dir, "output");
      try (FileOutputStream fos = new FileOutputStream(output)) {
        try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(sdf)) {
          NBlockDetector.detectNs(reader, 3, fos);
        }
      }
      assertEquals(EXPECTED, FileUtils.fileToString(output));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }


}
