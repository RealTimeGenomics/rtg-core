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
import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class AlternatingSequencesWriterTest extends TestCase {


  public void testAlternation() throws Exception {

    final File tempDir = FileUtils.createTempDir("cg2sdftest", "version2");
    try {

      final File reads = new File(tempDir, "reads");
      FileHelper.resourceToFile("com/rtg/reader/resources/sample.tsv", reads);

      final AlternatingSequencesWriter asw = new AlternatingSequencesWriter(new TsvSequenceDataSource(reads, 5), tempDir, 1000, PrereadType.CG, false, null);

      asw.processSequences(true, false);
      assertEquals(10, asw.mMaxLength);
      assertEquals(10, asw.mMinLength);
      assertEquals(6, asw.mNumberOfSequences);
      assertEquals(60, asw.mTotalLength);

      final MemoryPrintStream mps = new MemoryPrintStream();

      new SdfStatistics().mainInit(new String[] {tempDir.getPath() + "/left"}, mps.outputStream(), mps.printStream());
      TestUtils.containsAll(mps.toString(), "",
        "DNA",
        "CG",
        "LEFT",
        "Number of sequences: 3",
        "Maximum length     : 10",
        "Minimum length     : 10",
        "N                  : 0",
        "A                  : 15",
        "C                  : 0",
        "G                  : 10",
        "T                  : 5",
        "Total residues     : 30",
        "Residue qualities  : yes");


      mps.reset();

      new SdfStatistics().mainInit(new String[] {tempDir.getPath() + "/right"}, mps.outputStream(), mps.printStream());
      TestUtils.containsAll(mps.toString(), "",
        "DNA",
        "CG",
        "RIGHT",
        "Number of sequences: 3",
        "Maximum length     : 10",
        "Minimum length     : 10",
        "N                  : 0",
        "A                  : 5",
        "C                  : 10",
        "G                  : 0",
        "T                  : 15",
        "Total residues     : 30",
        "Residue qualities  : yes");

    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private static final String SEQ_L = "AAACCCGGGTTT";
  private static final String SEQ_R = "TTTGGGCCCAAA";
  private static final String QUAL_L = "555555555555";
  private static final String QUAL_R = "555555555555";
  private static final String TEST_TSV = ""
          + "0\t" + SEQ_L + SEQ_R + "\t" + QUAL_L + QUAL_R + StringUtils.LS;

  private void checkInMemory(final LongRange range) throws IOException {
    final File dir = FileUtils.createTempDir("test", "directtsv");
    try {
      final File tsvF = FileUtils.stringToFile(TEST_TSV, new File(dir, "cg.tsv"));
      final TsvSequenceDataSource tsv = new TsvSequenceDataSource(tsvF, 5);
      final AlternatingSequencesWriter alt = new AlternatingSequencesWriter(tsv, null, PrereadType.CG, true);
      final CompressedMemorySequencesReader[] rs = alt.processSequencesInMemoryPaired(tsvF, true, null, null, range);
      assertTrue(rs[0].nextSequence());
      final byte[] dna = new byte[(int) rs[0].maxLength()];
      rs[0].readCurrent(dna);
      assertEquals(SEQ_L, DnaUtils.bytesToSequenceIncCG(dna));
      assertTrue(rs[1].nextSequence());
      rs[1].readCurrent(dna);
      assertEquals(SEQ_R, DnaUtils.bytesToSequenceIncCG(dna));
      assertFalse(rs[0].nextSequence());
      assertFalse(rs[1].nextSequence());
      assertEquals(PrereadType.CG, rs[0].getPrereadType());
      assertEquals(PrereadType.CG, rs[1].getPrereadType());
      assertEquals(tsvF, rs[0].path());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testInMemory() throws IOException {
    checkInMemory(LongRange.NONE);
  }

  public void testInMemoryWithNoEndRange() throws IOException {
    checkInMemory(new LongRange(0, LongRange.MISSING));
  }

  public void testInMemoryWithRange() throws IOException {
    checkInMemory(new LongRange(0, 1));
  }
}
