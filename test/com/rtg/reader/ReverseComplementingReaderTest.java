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
import java.io.InputStream;
import java.util.ArrayList;

import com.rtg.mode.DNA;
import com.rtg.mode.DNAFastaSymbolTable;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class ReverseComplementingReaderTest extends DefaultSequencesReaderTest {

  public static Test suite() {
    return new TestSuite(ReverseComplementingReaderTest.class);
  }

  @Override
  protected SequencesReader createSequencesReader(final File dir) throws IOException {
    return new ReverseComplementingReader(SequencesReaderFactory.createDefaultSequencesReader(dir));
  }

  @Override
  protected byte[] getExpectedQuality() {
    final byte[] r = super.getExpectedQuality();
    ReverseComplementingReader.reverse(r, 0, r.length);
    return r;
  }

  // Can't do protein stuff with this reader

  @Override
  public void testEquals() { }

  @Override
  public void testResidueCountsWithoutDusterOnProtein() { }

  @Override
  public void testResidueCountsWithDusterOnProtein() { }

  @Override
  public void testReadRoll() { }

  @Override
  public void testRawRead() { }

  @Override
  public void testReadQualityMultifile() { }

  @Override
  public void testRead() { }

  @Override
  public void testResidueCountsWithDusterOnDNA() throws Exception {
    //create data source
    ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nac\n  tg\ntnGh\n\n\t   \n>test2\r\nATGC"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
    ds.setDusting(true);
    SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      final SequencesIterator it = dsr.iterator();
      assertEquals(mDir, dsr.path());
      assertTrue(it.nextSequence());
      assertEquals(0, it.currentSequenceId());
      SequencesWriterTest.checkEquals(it, new byte[]{0, 2, 0, 0, 0, 0, 0, 0});
      assertTrue(it.nextSequence());
      assertEquals(1, it.currentSequenceId());
      SequencesWriterTest.checkEquals(it, new byte[]{3, 2, 1, 4});
      assertEquals(1, dsr.residueCounts()[DNA.T.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.C.ordinal()]);
      assertEquals(1, dsr.residueCounts()[DNA.G.ordinal()]);
      assertEquals(1, dsr.residueCounts()[DNA.A.ordinal()]);
      assertEquals(7, dsr.residueCounts()[DNA.N.ordinal()]);
      assertFalse(it.nextSequence());
    }
  }

  @Override
  public void testResidueCountsWithoutDusterOnDNA() throws Exception {
    //create data source
    ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test\nac\n  tg\ntnGh\n\n\t   \n>test2\r\nATGC"));
    FastaSequenceDataSource ds = new FastaSequenceDataSource(al,  new DNAFastaSymbolTable());
    SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

    //testing the read (stolen from SequencesWriterTest)
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      final SequencesIterator it = dsr.iterator();
      assertEquals(mDir, dsr.path());
      assertEquals(12, dsr.totalLength());
      assertTrue(it.nextSequence());
      assertEquals(0, it.currentSequenceId());
      SequencesWriterTest.checkEquals(it, new byte[]{0, 2, 0, 1, 2, 1, 3, 4});
      assertTrue(it.nextSequence());
      assertEquals(1, it.currentSequenceId());
      SequencesWriterTest.checkEquals(it, new byte[]{3, 2, 1, 4});
      assertEquals(2, dsr.residueCounts()[DNA.N.ordinal()]);
      assertEquals(3, dsr.residueCounts()[DNA.A.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.T.ordinal()]);
      assertEquals(3, dsr.residueCounts()[DNA.C.ordinal()]);
      assertEquals(2, dsr.residueCounts()[DNA.G.ordinal()]);
      assertFalse(it.nextSequence());
    }
  }

  public void testRC() throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(createStream(">test1\nacngta\n"));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20000, PrereadType.UNKNOWN, false);
    sw.processSequences();
    try (SequencesReader dsr = createSequencesReader(mDir)) {
      final SequencesIterator it = dsr.iterator();
      assertTrue(it.nextSequence());
      final byte[] x = new byte[27];
      assertEquals(6, it.readCurrent(x));
      assertEquals(DNA.T.ordinal(), x[0]);
      assertEquals(DNA.A.ordinal(), x[1]);
      assertEquals(DNA.C.ordinal(), x[2]);
      assertEquals(DNA.N.ordinal(), x[3]);
      assertEquals(DNA.G.ordinal(), x[4]);
      assertEquals(DNA.T.ordinal(), x[5]);
      assertFalse(it.nextSequence());
    }
  }
}

