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
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Arrays;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test class
 */
public class CompressedMemorySequencesWriterTest extends TestCase {
  private static final String SEQUENCE_1 = "acgtacgtacgtacgt";
  private static final String QUALITY_1 =  "%^&*%^&*%^&*%^&*";
  private static final String SEQUENCE_2 = "cgatcgatcgatcgatcgat";
  private static final String QUALITY_2 =  "%^&*%^&*%^&*%^&*%^&*";
  private static final String FASTA = ">sequence_the_first\n" + SEQUENCE_1 + "\n>sequence_the_second\n" + SEQUENCE_2;
  private static final String FASTQ = "@sequence_the_first\n" + SEQUENCE_1 + "\n+sequence_the_first\n" + QUALITY_1 + "\n@sequence_the_second\n" + SEQUENCE_2 + "\n+sequence_the_second\n" + QUALITY_2;

  public void testFasta() throws IOException {
    SequencesWriter sw = new SequencesWriter(getFastaSource(FASTA), null, PrereadType.UNKNOWN, true);
    CompressedMemorySequencesReader cmsr = sw.processSequencesInMemory(null, true, new SimplePrereadNames(), new SimplePrereadNames(), LongRange.NONE);
    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_first", cmsr.currentName());
    assertEquals(16L, cmsr.currentLength());
    assertEquals(0L, cmsr.currentSequenceId());
    byte[] exp = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    byte[] res = new byte[16];
    cmsr.readCurrent(res);
    assertTrue("Exp: " + Arrays.toString(exp) + "\nact: " + Arrays.toString(res), Arrays.equals(exp, res));
    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_second", cmsr.currentName());
    assertEquals(20L, cmsr.currentLength());
    assertEquals(1L, cmsr.currentSequenceId());
    exp = new byte[] {2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4};
    res = new byte[20];
    cmsr.readCurrent(res);
    assertTrue("Exp: " + Arrays.toString(exp) + "\nact: " + Arrays.toString(res), Arrays.equals(exp, res));
    assertFalse(cmsr.nextSequence());
  }

  public void testFastq() throws IOException {
    SequencesWriter sw = new SequencesWriter(getFastqSource(FASTQ), null, PrereadType.UNKNOWN, true);
    CompressedMemorySequencesReader cmsr = sw.processSequencesInMemory(null, true, new SimplePrereadNames(), new SimplePrereadNames(), LongRange.NONE);
    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_first", cmsr.currentName());
    assertEquals(16L, cmsr.currentLength());
    assertEquals(0L, cmsr.currentSequenceId());
    byte[] exp = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    byte[] res = new byte[16];
    cmsr.readCurrent(res);
    assertTrue("Exp: " + Arrays.toString(exp) + "\nact: " + Arrays.toString(res), Arrays.equals(exp, res));
    byte[] expQual = {4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9};
    byte[] resQual = new byte[16];
    cmsr.readCurrentQuality(resQual);
    assertTrue("Exp: " + Arrays.toString(expQual) + "\nact: " + Arrays.toString(resQual), Arrays.equals(expQual, resQual));
    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_second", cmsr.currentName());
    assertEquals(20L, cmsr.currentLength());
    assertEquals(1L, cmsr.currentSequenceId());
    exp = new byte[] {2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4, 2, 3, 1, 4};
    res = new byte[20];
    cmsr.readCurrent(res);
    assertTrue("Exp: " + Arrays.toString(exp) + "\nact: " + Arrays.toString(res), Arrays.equals(exp, res));
    expQual = new byte[] {4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9};
    resQual = new byte[20];
    cmsr.readCurrentQuality(resQual);
    assertTrue("Exp: " + Arrays.toString(expQual) + "\nact: " + Arrays.toString(resQual), Arrays.equals(expQual, resQual));

    assertFalse(cmsr.nextSequence());
  }

  public void testFastqRegion() throws IOException {
    SequencesWriter sw = new SequencesWriter(getFastqSource(FASTQ), null, PrereadType.UNKNOWN, true);
    CompressedMemorySequencesReader cmsr = sw.processSequencesInMemory(null, true, new SimplePrereadNames(), new SimplePrereadNames(), new LongRange(0, 1));
    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_first", cmsr.currentName());
    assertEquals(16L, cmsr.currentLength());
    assertEquals(0L, cmsr.currentSequenceId());
    byte[] exp = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    byte[] res = new byte[16];
    cmsr.readCurrent(res);
    assertTrue("Exp: " + Arrays.toString(exp) + "\nact: " + Arrays.toString(res), Arrays.equals(exp, res));
    byte[] expQual = {4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9, 4, 61, 5, 9};
    byte[] resQual = new byte[16];
    cmsr.readCurrentQuality(resQual);
    assertTrue("Exp: " + Arrays.toString(expQual) + "\nact: " + Arrays.toString(resQual), Arrays.equals(expQual, resQual));
    assertFalse(cmsr.nextSequence());
  }

  public void testFastqRegionLarge() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(baos));
    SequencesWriter sw = new SequencesWriter(getFastqSource(FASTQ), null, PrereadType.UNKNOWN, true);
    CompressedMemorySequencesReader cmsr = sw.processSequencesInMemory(null, true, new SimplePrereadNames(), new SimplePrereadNames(), new LongRange(0, 10));

    assertTrue("missing out of range message", baos.toString().contains("The end sequence id \"10\" is out of range, it must be from \"1\" to \"2\". Defaulting end to \"2\""));

    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_first", cmsr.currentName());
    assertEquals(16L, cmsr.currentLength());
    assertEquals(0L, cmsr.currentSequenceId());

    assertTrue(cmsr.nextSequence());
    assertEquals("sequence_the_second", cmsr.currentName());
    assertEquals(20L, cmsr.currentLength());
    assertEquals(1L, cmsr.currentSequenceId());

    assertFalse(cmsr.nextSequence());
  }

  private FastaSequenceDataSource getFastaSource(String str) {
    return new FastaSequenceDataSource(Arrays.asList(new InputStream[] {new ByteArrayInputStream(str.getBytes())}), new DNAFastaSymbolTable());
  }


  private FastqSequenceDataSource getFastqSource(String str) {
    return new FastqSequenceDataSource(Arrays.asList(new InputStream[] {new ByteArrayInputStream(str.getBytes())}), FastqSequenceDataSource.FastQScoreType.PHRED);
  }
}
