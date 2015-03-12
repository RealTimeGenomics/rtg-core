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
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import com.rtg.util.StringUtils;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class CompressedMemorySequencesReader2Test extends CompressedMemorySequencesReaderTest {

  @Override
  protected SequencesReader createSequencesReader(final File dir, LongRange region) throws IOException {
    return CompressedMemorySequencesReader.createSequencesReader(dir, true, false, region, true);
  }

  private static final String SEQ_DATA = ""
      + "@read0" + StringUtils.LS
      + "GCGTGAGTACGTGACTGAGCGGCATGCTGAAATCC" + StringUtils.LS           //35
      + "+read0" + StringUtils.LS
      + "OSQSSRNS'&*DSQNPM%KPO$PQRHHHRNPQIRS" + StringUtils.LS
      + "@read1" + StringUtils.LS
      + "AAATCG" + StringUtils.LS                                        //6
      + "+read1" + StringUtils.LS
      + "ERQSRI" + StringUtils.LS
      + "@read2" + StringUtils.LS
      + "GAGGATCCTTAAGTGTC" + StringUtils.LS                             //17
      + "+read2" + StringUtils.LS
      + "OHPPGRQSORFRSRQSP" + StringUtils.LS
      + "@read3" + StringUtils.LS
      + "ACTTTTCAGCAGGTCGAGTTCTCCTATACTGAAAG" + StringUtils.LS           //35
      + "+read3" + StringUtils.LS
      + "SQNQNMPNIPLIPGQRI?PMPOOSMRISSRHPMSR" + StringUtils.LS
      + "@read4" + StringUtils.LS
      + "ACTCGTCTAGTACGTTGTCCACGAACCAGG" + StringUtils.LS                //30
      + "+read4" + StringUtils.LS
      + "SLOODPNLR=$NNRQOPSRINRNRSPSOMR" + StringUtils.LS;

  public void testData() throws IOException {
    final File dir = FileUtils.createTempDir("test", "dir");
    try {
      ReaderTestUtils.getReaderDNAFastq(SEQ_DATA, dir, new SdfId(), 20);
      final IndexFile f = new IndexFile(dir);

      final CompressedMemorySequencesReader2 mem = new CompressedMemorySequencesReader2(dir, f, true, true, new LongRange(0, 5));
      final FastqSequenceDataSource fastq = new FastqSequenceDataSource(Arrays.asList((InputStream) new ByteArrayInputStream(SEQ_DATA.getBytes())), FastqSequenceDataSource.FastQScoreType.PHRED);
      int i = 0;
      while (fastq.nextSequence()) {
        final byte[] exp = Arrays.copyOf(fastq.sequenceData(), fastq.currentLength());
        final byte[] act = new byte[mem.length(i)];
        assertEquals(exp.length, act.length);
        mem.read(i, act);
        assertTrue("i: " + i + " \nExp: " + Arrays.toString(exp) + " \nACT: " + Arrays.toString(act), Arrays.equals(exp, act));
        final byte[] expQual = Arrays.copyOf(fastq.qualityData(), fastq.currentLength());
        final byte[] actQual = new byte[act.length];
        mem.readQuality(i, actQual);
        assertTrue(Arrays.equals(expQual, actQual));

        i++;
      }
      final SequencesReader mem2 = mem.copy();
      final FastqSequenceDataSource fastq2 = new FastqSequenceDataSource(Arrays.asList((InputStream) new ByteArrayInputStream(SEQ_DATA.getBytes())), FastqSequenceDataSource.FastQScoreType.PHRED);
      int i2 = 0;
      while (fastq2.nextSequence()) {
        final byte[] exp = Arrays.copyOf(fastq2.sequenceData(), fastq2.currentLength());
        final byte[] act = new byte[mem2.length(i2)];
        assertEquals(exp.length, act.length);
        mem.read(i2, act);
        assertTrue("i: " + i2 + " \nExp: " + Arrays.toString(exp) + " \nACT: " + Arrays.toString(act), Arrays.equals(exp, act));
        final byte[] expQual = Arrays.copyOf(fastq2.qualityData(), fastq2.currentLength());
        final byte[] actQual = new byte[act.length];
        mem2.readQuality(i2, actQual);
        assertTrue(Arrays.equals(expQual, actQual));
        i2++;
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testLastNameBug() throws IOException {
    final File dir = FileUtils.createTempDir("test", "dir");
    try {
      final String data = FileHelper.resourceToString("com/rtg/reader/resources/reads100.fastq");
      ReaderTestUtils.getReaderDNAFastq(data, dir, new SdfId(), 1234);
      final IndexFile f = new IndexFile(dir);

      final LongRange region = new LongRange(0, 1);
      final CompressedMemorySequencesReader2 mem = new CompressedMemorySequencesReader2(dir, f, true, true, region);
      final CompressedMemorySequencesReader memold = new CompressedMemorySequencesReader(dir, f, 5, true, true, region);

      for (int i = 0; i < memold.numberSequences(); i++) {
        final String oldname = memold.fullName(i);
        final String name = mem.fullName(i);
        assertEquals("i: " + i, oldname, name);
      }
      try {
        mem.fullName(memold.numberSequences());
        fail("Got name beyond number of sequences.");
      } catch (IndexOutOfBoundsException e) {
      }

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  //
//  public void testFoo() {
//     final File f = new File("/rtgshare/data/human/ref/1000g_v37_phase1/sdf");
//    SequencesReader c = CompressedMemorySequencesReader.createSequencesReader(f, true, false, LongRange.NONE);
//    byte[] buf = new byte[10000000];
//        c.seek(4);
//        int pos = 120000000;
//        while (pos < c.currentLength()) {
//          try {
//          pos += c.readCurrent(buf, pos, Math.min(buf.length, c.currentLength() - pos));
//          } catch (RuntimeException e) {
//            System.err.println(c.currentName() + " " + pos);
//            throw e;
//          }
//        }
//  }

}
