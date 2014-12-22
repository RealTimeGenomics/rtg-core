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
        final SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(sdf);
        try {
          NBlockDetector.detectNs(reader, 3, fos);
        } finally {
          reader.close();
        }
      }
      assertEquals(EXPECTED, FileUtils.fileToString(output));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }


}
