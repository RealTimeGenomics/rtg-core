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
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.Protein;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.reader.FastqSequenceDataSource.FastQScoreType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for <code>SequencesReader</code> implementations
 */
public abstract class AbstractSequencesReaderTest extends TestCase {

  protected abstract SequencesReader createSequencesReader(final File dir) throws IOException;
  protected abstract boolean canReadTwice();

  public AbstractSequencesReaderTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(AbstractSequencesReaderTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  protected File mDir;

  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws Exception {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  protected InputStream createStream(final String data) {
    return new ByteArrayInputStream(data.getBytes());
  }

  public void testLabel() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgt\n>bob\ntagt\n>hobos r us\naccc"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 30, PrereadType.UNKNOWN, false);
    sw.processSequences();
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(Exam.integrity(dsr));

      assertEquals(mDir, dsr.path());
      try {
        dsr.currentName();
        fail();
      } catch (final IllegalStateException e) {
        assertEquals("Last call to nextSequence() or seek() failed and left current information unavailable.", e.getMessage());
      }
      assertTrue(dsr.nextSequence());
      assertEquals("test", dsr.currentName());
      assertEquals("test", dsr.currentName());
      assertTrue(dsr.nextSequence());
      assertEquals("bob", dsr.currentName());
      assertTrue(dsr.nextSequence());
      assertEquals("hobos", dsr.currentName());
      dsr.seek(1);
      assertEquals("bob", dsr.currentName());
      dsr.seek(0);
      assertEquals("test", dsr.currentName());
      dsr.seek(2);
      assertEquals("hobos", dsr.currentName());
    }
  }

  public void testLengthBetween() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    //32
    //17
    //20
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(32, dsr.lengthBetween(0, 1));
      assertEquals(17, dsr.lengthBetween(1, 2));
      assertEquals(20, dsr.lengthBetween(2, 3));
      assertEquals(49, dsr.lengthBetween(0, 2));
      assertEquals(37, dsr.lengthBetween(1, 3));
      assertEquals(69, dsr.lengthBetween(0, 3));
      assertEquals(0, dsr.lengthBetween(3, 3));
    }
  }

  public void testLengthBetween2() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nactggtcatgca\n>bob the buuilder\ntagttcagcatc\n>hobos r us\naccccaccccac"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(12, dsr.lengthBetween(0, 1));
      assertEquals(12, dsr.lengthBetween(1, 2));
      assertEquals(12, dsr.lengthBetween(2, 3));
      assertEquals(24, dsr.lengthBetween(0, 2));
      assertEquals(24, dsr.lengthBetween(1, 3));
      assertEquals(36, dsr.lengthBetween(0, 3));
      assertEquals(0, dsr.lengthBetween(3, 3));
    }
  }

  public void testLengths() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    //32
    //17
    //20
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(Arrays.equals(new int[]{32, 17, 20}, dsr.sequenceLengths(0, 3)));
      assertTrue(Arrays.equals(new int[]{17, 20}, dsr.sequenceLengths(1, 3)));
      assertTrue(Arrays.equals(new int[]{20}, dsr.sequenceLengths(2, 3)));
      assertTrue(Arrays.equals(new int[]{32, 17}, dsr.sequenceLengths(0, 2)));
      assertTrue(Arrays.equals(new int[]{17}, dsr.sequenceLengths(1, 2)));
      assertTrue(Arrays.equals(new int[]{32}, dsr.sequenceLengths(0, 1)));
    }
  }

  public void testLengths2() throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    //32
    //17
    //20
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 2000000, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(Arrays.equals(new int[]{32, 17, 20}, dsr.sequenceLengths(0, 3)));
      assertTrue(Arrays.equals(new int[]{17, 20}, dsr.sequenceLengths(1, 3)));
      assertTrue(Arrays.equals(new int[]{20}, dsr.sequenceLengths(2, 3)));
      assertTrue(Arrays.equals(new int[]{32, 17}, dsr.sequenceLengths(0, 2)));
      assertTrue(Arrays.equals(new int[]{17}, dsr.sequenceLengths(1, 2)));
      assertTrue(Arrays.equals(new int[]{32}, dsr.sequenceLengths(0, 1)));
    }
  }

  public void testReadRoll() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();


    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(Exam.integrity(dsr));

      assertEquals(mDir, dsr.path());
      try {
        dsr.currentName();
        fail();
      } catch (final IllegalStateException e) {
        assertEquals("Last call to nextSequence() or seek() failed and left current information unavailable.", e.getMessage());
      }
      assertTrue(dsr.nextSequence());
      assertEquals(0, dsr.currentSequenceId());
      assertEquals("test", dsr.currentName());
      assertEquals(32, dsr.currentLength());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 2, 3, 4, 3, 4, 3, 4, 3, 4, 2, 4, 4, 1, 3, 3, 3, 2, 4, 2, 1, 2, 4, 3, 3, 4, 2, 1, 4, 3, 2, 1});
      assertTrue(dsr.nextSequence());
      assertEquals(1, dsr.currentSequenceId());
      assertEquals("bob", dsr.currentName());
      assertEquals(17, dsr.currentLength());
      SequencesWriterTest.checkEquals(dsr, new byte[]{4, 1, 3, 4, 4, 2, 1, 3, 2, 1, 4, 2, 3, 1, 4, 2, 1});
      assertTrue(Exam.integrity(dsr));
      assertTrue(dsr.nextSequence());
      assertEquals(2, dsr.currentSequenceId());
      assertEquals("hobos", dsr.currentName());
      assertEquals(20, dsr.currentLength());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1});
      dsr.seek(1);
      assertEquals(1, dsr.currentSequenceId());
      assertEquals("bob", dsr.currentName());
      assertEquals(17, dsr.currentLength());
      SequencesWriterTest.checkEquals(dsr, new byte[]{4, 1, 3, 4, 4, 2, 1, 3, 2, 1, 4, 2, 3, 1, 4, 2, 1});
      dsr.reset();
      final byte[] seqs = new byte[32];
      final int amount = dsr.read(0, seqs);
      DNA[] expected = {DNA.A, DNA.C, DNA.G, DNA.T, DNA.G, DNA.T, DNA.G, DNA.T,
        DNA.G, DNA.T, DNA.C, DNA.T, DNA.T, DNA.A, DNA.G, DNA.G, DNA.G, DNA.C, DNA.T,
        DNA.C, DNA.A, DNA.C, DNA.T, DNA.G, DNA.G, DNA.T, DNA.C, DNA.A, DNA.T, DNA.G,
        DNA.C, DNA.A};
      assertEquals(expected.length, amount);
      for (int i = 0; i < seqs.length; i++) {
        assertEquals(expected[i].ordinal(), seqs[i]);
      }
      assertTrue(Exam.integrity(dsr));
      dsr.seek(0);
      expected = new DNA[]{DNA.A, DNA.C, DNA.G, DNA.T, DNA.G, DNA.T, DNA.G, DNA.T,
        DNA.G, DNA.T, DNA.C, DNA.T, DNA.T, DNA.A, DNA.G, DNA.G, DNA.G, DNA.C, DNA.T,
        DNA.C, DNA.A, DNA.C, DNA.T, DNA.G, DNA.G, DNA.T, DNA.C, DNA.A, DNA.T, DNA.G,
        DNA.C, DNA.A};
      byte[] bytes = new byte[expected.length];
      dsr.readCurrent(bytes);
      DNA[] dnaValues = DNA.values();
      for (int i = 0; i < bytes.length; i++) {
        assertEquals(expected[i], dnaValues[bytes[i]]);
      }

      dsr.seek(0);
      expected = new DNA[]{DNA.G, DNA.T, DNA.G, DNA.T,
        DNA.G, DNA.T, DNA.C, DNA.T, DNA.T, DNA.A, DNA.G, DNA.G, DNA.G, DNA.C, DNA.T,
        DNA.C, DNA.A, DNA.C, DNA.T, DNA.G, DNA.G, DNA.T, DNA.C, DNA.A, DNA.T, DNA.G,
        DNA.C, DNA.A};
      bytes = new byte[expected.length];
      dsr.readCurrent(bytes, 4, expected.length);
      dnaValues = DNA.values();
      for (int i = 0; i < bytes.length; i++) {
        assertEquals(expected[i], dnaValues[bytes[i]]);
      }

    }
  }

  public void testRead() throws Exception {
    //create data source
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nac\n  tg\ntnGh\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    //testing the read (stolen from SequencesWriterTest)
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      assertTrue(dsr.nextSequence());
      assertEquals(0, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 2, 4, 3, 4, 0, 3, 0});
      assertTrue(dsr.nextSequence());
      assertEquals(1, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 4, 3, 2});
      assertFalse(dsr.nextSequence());
    }
  }

  protected byte[] getExpectedQuality() {
    return new byte[] {'!' - '!', '<' - '!', '>' - '!', '<' - '!', '#' - '!', '#' - '!', '!' - '!', '<' - '!'};
  }

  public void testReadQuality() throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream("@testQuality\n"
                        + "actgcatc\n"
                        + "+\n"
                        + "!<><##!<"));
    final FastqSequenceDataSource fq = new FastqSequenceDataSource(al, FastQScoreType.PHRED);
    final SequencesWriter sw = new SequencesWriter(fq, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(dsr.nextSequence());
      assertEquals("testQuality", dsr.currentName());
      final byte[] qual = new byte[dsr.currentLength()];
      assertEquals(qual.length, dsr.readCurrentQuality(qual));
      final byte[] exp = getExpectedQuality();

      assertTrue(Arrays.equals(exp, qual));

      dsr.reset();
      assertEquals(qual.length, dsr.readQuality(0, qual));
      assertTrue(Arrays.equals(exp, qual));

      if (!canReadTwice()) {
        try {
          dsr.readCurrentQuality(qual);
          fail();
        } catch (final IllegalStateException e) {
          // correct
        }
      }
      dsr.seek(0);
      assertFalse(dsr.nextSequence());
      try {
        dsr.readCurrentQuality(qual);
        fail();
      } catch (final IllegalStateException e) {
        // correct
      }
    }
  }

  public void testReadQualityMultifile() throws IOException {
    //check rolling works
    final byte[] expQualBytes = new byte[60];
    final StringBuilder quals = new StringBuilder();
    for (int i = 0; i < expQualBytes.length; i++) {
      expQualBytes[i] = (byte) (60 - i);
      quals.append((char) (60 - i + '!'));
    }
    final byte[] expDnaBytes = new byte[60];
    final StringBuilder dnas = new StringBuilder();
    for (int i = 0; i < expDnaBytes.length; i++) {
      expDnaBytes[i] = (byte) DNA.A.ordinal();
      dnas.append(DNA.A.toString());
    }
    final String seq = "@TotallyUniqueName\n" + dnas.toString() + "\n+\n" + quals.toString();
    final ArrayList<InputStream> streams = new ArrayList<>();
    streams.add(createStream(seq));
    final FastqSequenceDataSource ds = new FastqSequenceDataSource(streams, FastQScoreType.PHRED);
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertTrue(dsr.nextSequence());
      final byte[] result = new byte[60];
      assertEquals(60, dsr.currentLength());
      assertEquals("TotallyUniqueName", dsr.currentName());
      assertEquals(60, dsr.readCurrentQuality(result));
      assertTrue(Arrays.equals(expQualBytes, result));
    }
  }

  protected DNA[] getExpected0() {
    return new DNA[] {DNA.A, DNA.C, DNA.G, DNA.T, DNA.G, DNA.T, DNA.G, DNA.T,
                                   DNA.G, DNA.T, DNA.C, DNA.T, DNA.T, DNA.A, DNA.G, DNA.G, DNA.G, DNA.C, DNA.T,
                                   DNA.C, DNA.A, DNA.C, DNA.T, DNA.G, DNA.G, DNA.T, DNA.C, DNA.A, DNA.T, DNA.G,
                                   DNA.C, DNA.A};
  }

  protected DNA[] getExpected1() {
    return new DNA[] {DNA.T, DNA.A, DNA.G, DNA.T, DNA.T, DNA.C, DNA.A, DNA.G, DNA.C,
                                   DNA.A, DNA.T, DNA.C, DNA.G, DNA.A, DNA.T, DNA.C, DNA.A};
  }

  protected DNA[] getExpected2() {
    return new DNA[] {DNA.A, DNA.C, DNA.C, DNA.C, DNA.C, DNA.A, DNA.C, DNA.C, DNA.C,
                      DNA.C, DNA.A, DNA.C, DNA.A, DNA.A, DNA.A, DNA.C, DNA.C, DNA.C, DNA.A, DNA.A};
  }

  public void testRawRead() throws Exception {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();


    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      final byte[] bigenough = new byte[32];
      final byte[] notbigenough = new byte[19];
      final DNA[] expected0 = getExpected0();
      final DNA[] expected1 = getExpected1();
      final DNA[] expected2 = getExpected2();
      dsr.read(0, bigenough);
      for (int i = 0; i < expected0.length; i++) {
        assertEquals(expected0[i].ordinal(), bigenough[i]);
      }
      dsr.read(1, bigenough);
      for (int i = 0; i < expected1.length; i++) {
        assertEquals(expected1[i].ordinal(), bigenough[i]);
      }
      dsr.read(1, bigenough, 8, expected1.length - 8);
      for (int i = 8; i < expected1.length; i++) {
        assertEquals(String.valueOf(i), expected1[i].ordinal(), bigenough[i - 8]);
      }
      dsr.read(2, bigenough);
      for (int i = 0; i < expected2.length; i++) {
        assertEquals(expected2[i].ordinal(), bigenough[i]);
      }
      dsr.read(1, notbigenough);
      for (int i = 0; i < expected1.length; i++) {
        assertEquals(expected1[i].ordinal(), notbigenough[i]);
      }
      try {
        dsr.read(2, notbigenough);
        fail("Should have complained about array size");
      } catch (final IllegalArgumentException e) {
        assertEquals("Array too small got: " + notbigenough.length + " required: 20", e.getMessage());
      }
      try {
        dsr.read(0, notbigenough);
        fail("Should have complained about array size");
      } catch (final IllegalArgumentException e) {
        assertEquals("Array too small got: " + notbigenough.length + " required: 32", e.getMessage());
      }
    }
  }

  public void testRead2() throws Exception {
    final ArrayList<InputStream> al1 = new ArrayList<>();
    al1.add(createStream(">x1\nACGTACGNNNNN\n"));
    final FastaSequenceDataSource ds1 = new FastaSequenceDataSource(al1,
                                                              new DNAFastaSymbolTable());
    final SequencesWriter sw1 = new SequencesWriter(ds1, mDir, 20, PrereadType.UNKNOWN, false);
    sw1.processSequences();

    try (SequencesReader dsr2 = createSequencesReader(mDir)) {
      assertEquals(dsr2.totalLength(), 12);

      dsr2.nextSequence();
      assertEquals(dsr2.totalLength(), dsr2.currentLength());
      assertEquals(dsr2.totalLength(), dsr2.maxLength());
      assertEquals(dsr2.totalLength(), dsr2.minLength());
      final byte[] bytes = new byte[(int) dsr2.maxLength()];
      dsr2.readCurrent(bytes);
      if (!canReadTwice()) {
        try {
          dsr2.readCurrent(bytes);
          fail("Should throw exception");
        } catch (final IllegalStateException e) {
          //
        }
      }
      assertTrue(!dsr2.nextSequence());
      try {
        dsr2.readCurrent(bytes);
        fail("Should throw exception.");
      } catch (final IllegalStateException e) {
        //
      }
    }
  }

  public void testEquals() throws Exception {
    final File dir1 = FileUtils.createTempDir("dir1", "11");
    final File dir2 = FileUtils.createTempDir("dir2", "22");
    //create data source
    final ArrayList<InputStream> al1 = new ArrayList<>();
    al1.add(createStream(">test\naH\n  tg\ntXGj\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds1 = new FastaSequenceDataSource(al1,
                                                              new ProteinFastaSymbolTable());
    final SequencesWriter sw1 = new SequencesWriter(ds1, dir1, 20, PrereadType.UNKNOWN, false);
    sw1.processSequences();


    final ArrayList<InputStream> al2 = new ArrayList<>();
    al2.add(createStream(">test\nacgtgtgtgtcttagggctcactggtcatgca\n>bob the buuilder\ntagttcagcatcgatca\n>hobos r us\naccccaccccacaaacccaa"));
    final FastaSequenceDataSource ds2 = new FastaSequenceDataSource(al2,
                                                              new DNAFastaSymbolTable());
    final SequencesWriter sw2 = new SequencesWriter(ds2, dir2, 20, PrereadType.UNKNOWN, false);
    sw2.processSequences();

    try (SequencesReader dsr1 = createSequencesReader(dir1)) {
      final SequencesReader dsr12 = createSequencesReader(dir1);
      try {
        final SequencesReader dsr2 = createSequencesReader(dir2);
        try {
          final SequencesReader dsr22 = createSequencesReader(dir2);
          try {
            assertTrue(!dsr1.equals(null));
            assertTrue(!dsr1.equals(dsr2));
            assertTrue(!dsr2.equals(dsr1));
            assertTrue(dsr2.equals(dsr22));
            assertEquals(dsr2.hashCode(), dsr22.hashCode());
          } finally {
            dsr22.close();
          }
        } finally {
          dsr2.close();
        }
      } finally {
        dsr12.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir1));
      assertTrue(FileHelper.deleteAll(dir2));
    }
  }

  public void testResidueCountsWithoutDusterOnDNA() throws Exception {
    //create data source
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nac\n  tg\ntnGh\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    //testing the read (stolen from SequencesWriterTest)
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      assertEquals(12, dsr.totalLength());
      assertTrue(dsr.nextSequence());
      assertEquals(0, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 2, 4, 3, 4, 0, 3, 0});
      assertTrue(dsr.nextSequence());
      assertEquals(1, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 4, 3, 2});
      assertEquals(2, dsr.residueCounts()[DNA.A.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.C.ordinal()]);
      assertEquals(3, dsr.residueCounts()[DNA.G.ordinal()]);
      assertEquals(3, dsr.residueCounts()[DNA.T.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.N.ordinal()]);
      assertFalse(dsr.nextSequence());
    }
  }

  public void testResidueCountsWithDusterOnDNA() throws Exception {
    //create data source
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nac\n  tg\ntnGh\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new DNAFastaSymbolTable());
    ds.setDusting(true);
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    //testing the read (stolen from SequencesWriterTest)
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      assertTrue(dsr.nextSequence());
      assertEquals(0, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{0, 0, 0, 0, 0, 0, 3, 0});
      assertTrue(dsr.nextSequence());
      assertEquals(1, dsr.currentSequenceId());
      SequencesWriterTest.checkEquals(dsr, new byte[]{1, 4, 3, 2});
      assertEquals(1, dsr.residueCounts()[DNA.A.ordinal()]);
      assertEquals(1, dsr.residueCounts()[DNA.C.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.G.ordinal()]);
      assertEquals(1, dsr.residueCounts()[DNA.T.ordinal()]);
      assertEquals(7, dsr.residueCounts()[DNA.N.ordinal()]);
      assertFalse(dsr.nextSequence());
    }
  }

  public void testResidueCountsWithoutDusterOnProtein() throws Exception {
    //create data source
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\naH\n  tg\ntXGj\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new ProteinFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    //testing the read (stolen from SequencesWriterTest)
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      assertTrue(dsr.nextSequence());
      SequencesWriterTest.checkEquals(dsr, new byte[]{2, 10, 18, 9, 18, 0, 9, 0});
      assertTrue(dsr.nextSequence());
      SequencesWriterTest.checkEquals(dsr, new byte[]{2, 18, 9, 6});
      assertEquals(2, dsr.residueCounts()[Protein.A.ordinal()]);
      assertEquals(1, dsr.residueCounts()[Protein.C.ordinal()]);
      assertEquals(1, dsr.residueCounts()[Protein.H.ordinal()]);
      assertEquals(3, dsr.residueCounts()[Protein.T.ordinal()]);
      assertEquals(3, dsr.residueCounts()[Protein.G.ordinal()]);
      assertEquals(2, dsr.residueCounts()[Protein.X.ordinal()]);
      assertFalse(dsr.nextSequence());
    }
  }

  public void testResidueCountsWithDusterOnProtein() throws Exception {
    //create data source
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\naH\n  tg\ntXGj\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
                                                             new ProteinFastaSymbolTable());
    ds.setDusting(true);
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    try (SequencesReader dsr = createSequencesReader(mDir)) {
      assertEquals(mDir, dsr.path());
      assertTrue(dsr.nextSequence());
      SequencesWriterTest.checkEquals(dsr, new byte[]{0, 10, 0, 0, 0, 0, 9, 0});
      assertTrue(dsr.nextSequence());
      SequencesWriterTest.checkEquals(dsr, new byte[]{2, 18, 9, 6});
      assertEquals(1, dsr.residueCounts()[Protein.A.ordinal()]);
      assertEquals(1, dsr.residueCounts()[Protein.C.ordinal()]);
      assertEquals(1, dsr.residueCounts()[Protein.H.ordinal()]);
      assertEquals(1, dsr.residueCounts()[Protein.T.ordinal()]);
      assertEquals(2, dsr.residueCounts()[Protein.G.ordinal()]);
      assertEquals(6, dsr.residueCounts()[Protein.X.ordinal()]);
      assertFalse(dsr.nextSequence());
    }
  }
}

