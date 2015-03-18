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
import com.rtg.util.TestUtils;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test for SdfFilter
 *
 */
public class SdfFilterTest extends TestCase {

  public SdfFilterTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(SdfFilterTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.runAndWait(suite());
  }

  private static byte[] makeDNA(final String s) {
    final byte[] dna = new byte[s.length()];
    for (int i = 0; i < dna.length; i++) {
      switch (s.charAt(i)) {
        case 'A':
        case 'a':
          dna[i] = (byte) DNA.A.ordinal();
          break;
        case 'C':
        case 'c':
          dna[i] = (byte) DNA.C.ordinal();
          break;
        case 'G':
        case 'g':
          dna[i] = (byte) DNA.G.ordinal();
          break;
        case 'T':
        case 't':
          dna[i] = (byte) DNA.T.ordinal();
          break;
        default:
          dna[i] = (byte) DNA.N.ordinal();
      }
    }
    return dna;
  }
  private InputStream createStream(final String data) {
    return new ByteArrayInputStream(data.getBytes());
  }

  public void testOneFile() throws IOException {
    final String testInputString = ">test without Ns\n"
        + "acgtgtgtgtcttagggctcactggtcatgca\n"
        + ">one N\n"
        + "tanttcagcatcgatca\n"
        + ">two Ns\n"
        + "ancccaccccanaaacccaa";

    final byte[][] expectedSequences = {
      //no Ns
      makeDNA("acgtgtgtgtcttagggctcactggtcatgca"),
      //1 N
      makeDNA("taattcagcatcgatca"),
      makeDNA("tacttcagcatcgatca"),
      makeDNA("tagttcagcatcgatca"),
      makeDNA("tatttcagcatcgatca"),
      //2 Ns
      makeDNA("aacccaccccanaaacccaa"),
      makeDNA("accccaccccanaaacccaa"),
      makeDNA("agcccaccccanaaacccaa"),
      makeDNA("atcccaccccanaaacccaa"),
    };
    final String[] expectedSequenceNames = {
      "test", "one-expN->A", "one-expN->C", "one-expN->G", "one-expN->T", "two-expN->A", "two-expN->C", "two-expN->G", "two-expN->T"
    };

    final File tempDir = FileUtils.createTempDir("sdffilter", "test");
    try {
      final File inDir = new File(tempDir, "test.in");
      final File outDir = new File(tempDir, "testout");
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(testInputString));
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[] {"-i", inDir.getPath(), "-o", outDir.getPath(), "-n", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());
      checkExpected(outDir, expectedSequences, expectedSequenceNames);
    } finally {
      assertTrue("Unable to delete " + tempDir + ": ", FileHelper.deleteAll(tempDir));
    }
  }
  public void testOneFileNoExpand() throws IOException {
    final String testInputString = ">test without Ns\n"
      + "acgtgtgtgtcttagggctcactggtcatgca\n"
      + ">one N\n"
      + "tanttcagcatcgatca\n"
      + ">two Ns\n"
      + "ancccaccccanaaacccaa";

    final byte[][] expectedSequences = {
        //no Ns
        makeDNA("acgtgtgtgtcttagggctcactggtcatgca"),
        //1 N
        makeDNA("tanttcagcatcgatca"),
        //2 Ns
        makeDNA("ancccaccccanaaacccaa"),
    };
    final String[] expectedSequenceNames = {
        "test", "one", "two"
    };

    final File tempDir = FileUtils.createTempDir("sdffilter", "test");
    try {
      final File inDir = new File(tempDir, "test.in");
      final File outDir = new File(tempDir, "testout");
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(testInputString));
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[]{"-i", inDir.getPath(), "-o", outDir.getPath()}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());
      checkExpected(outDir, expectedSequences, expectedSequenceNames);
    } finally {
      assertTrue("Unable to delete " + tempDir + ": ", FileHelper.deleteAll(tempDir));
    }
  }

  public void testPair() throws IOException {
    final String testInputString = ">test without Ns\n"
        + "acgtgtgtgtcttagggctcactggtcatgca\n"
        + ">one N\n"
        + "tanttcagcatcgatca\n"
        + ">two Ns\n"
        + "ancccaccccanaaacccaa";

    final byte[][] expectedSequences = {
      //no Ns
      makeDNA("acgtgtgtgtcttagggctcactggtcatgca"),
      //1 N
      makeDNA("taattcagcatcgatca"),
      makeDNA("tacttcagcatcgatca"),
      makeDNA("tagttcagcatcgatca"),
      makeDNA("tatttcagcatcgatca"),
      //2 Ns
      makeDNA("aacccaccccanaaacccaa"),
      makeDNA("accccaccccanaaacccaa"),
      makeDNA("agcccaccccanaaacccaa"),
      makeDNA("atcccaccccanaaacccaa"),
    };
    final String[] expectedSequenceNames = {
      "test", "one-expN->A", "one-expN->C", "one-expN->G", "one-expN->T", "two-expN->A", "two-expN->C", "two-expN->G", "two-expN->T"
    };

    final File dir = FileUtils.createTempDir("sdffilter", "testpair");
    try {
      final File leftInDir = new File(dir, "left-in");
      final File rightInDir = new File(dir, "right-in");
      final File outDir = new File(dir, "out");
      assertTrue(leftInDir.mkdir());
      assertTrue(rightInDir.mkdir());

      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(testInputString));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, leftInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();

      al = new ArrayList<>();
      al.add(createStream(testInputString));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, rightInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[] {"-l", leftInDir.getPath(), "-r", rightInDir.getPath(), "-o", outDir.getPath(), "-n", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());
      checkExpected(new File(outDir, "left"), expectedSequences, expectedSequenceNames);
    } finally {
      assertTrue("Unable to delete " + dir + ": ", FileHelper.deleteAll(dir));
    }
  }

  private void checkExpected(File outDir, byte[][] expectedSequences, String[] expectedSequenceNames) throws IOException {
    try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(outDir)) {
      assertEquals(expectedSequences.length, dsr.numberSequences());
      for (int j = 0; j < expectedSequences.length; j++) {
        assertEquals("j: " + j, expectedSequenceNames[j], dsr.name(j));
        assertTrue("j: " + j, Arrays.equals(expectedSequences[j], dsr.read(j)));
      }
    }
  }

  public void testPairAlternate() throws IOException {
    final String testInputString = ">test without Ns\n"
        + "acgtgtgtgtcttagggctcactggtcatgca\n"
        + ">one N\n"
        + "tanttcagcatcgatca\n"
        + ">two Ns\n"
        + "ancccaccccanaaacccaa";

    final byte[][] expectedSequences = {
      //no Ns
      makeDNA("acgtgtgtgtcttagggctcactggtcatgca"),
      //1 N left
      makeDNA("taattcagcatcgatca"),
      makeDNA("tacttcagcatcgatca"),
      makeDNA("tagttcagcatcgatca"),
      makeDNA("tatttcagcatcgatca"),
      //1 N right
      makeDNA("tanttcagcatcgatca"),
      makeDNA("tanttcagcatcgatca"),
      makeDNA("tanttcagcatcgatca"),
      makeDNA("tanttcagcatcgatca"),
      //2 Ns left
      makeDNA("ancccaccccaaaaacccaa"),
      makeDNA("ancccaccccacaaacccaa"),
      makeDNA("ancccaccccagaaacccaa"),
      makeDNA("ancccaccccataaacccaa"),
      //2 Ns right
      makeDNA("ancccaccccanaaacccaa"),
      makeDNA("ancccaccccanaaacccaa"),
      makeDNA("ancccaccccanaaacccaa"),
      makeDNA("ancccaccccanaaacccaa"),
    };
    final String[] expectedSequenceNames = {
      "test",
      "one-exp-left-N->A", "one-exp-left-N->C", "one-exp-left-N->G", "one-exp-left-N->T",
      "one-exp-right-N->A", "one-exp-right-N->C", "one-exp-right-N->G", "one-exp-right-N->T",
      "two-exp-left-N->A", "two-exp-left-N->C", "two-exp-left-N->G", "two-exp-left-N->T",
      "two-exp-right-N->A", "two-exp-right-N->C", "two-exp-right-N->G", "two-exp-right-N->T",
    };

    final File dir = FileUtils.createTempDir("sdffilter", "testpair");
    try {
      final File leftInDir = new File(dir, "left-in");
      final File rightInDir = new File(dir, "right-in");
      final File outDir = new File(dir, "out");
      assertTrue(leftInDir.mkdir());
      assertTrue(rightInDir.mkdir());

      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(testInputString));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, leftInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();

      al = new ArrayList<>();
      al.add(createStream(testInputString));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, rightInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[] {"-a", "-l", leftInDir.getPath(), "-r", rightInDir.getPath(), "-o", outDir.getPath(), "-n", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());

      checkExpected(new File(outDir, "left"), expectedSequences, expectedSequenceNames);
    } finally {
      assertTrue("Unable to delete " + dir + ": ", FileHelper.deleteAll(dir));
    }
  }

  public void testPairAlternate2() throws IOException {
    final String withN = "acgtacntacgt";
    final String testInputStringWithN = ">test\n" + withN;
    final String withoutN = "acgtacgtacgt";
    final String testInputStringWithoutN = ">test\n" + withoutN;
    final byte[][] expectedSequencesWithN = {
      makeDNA("acgtacatacgt"),
      makeDNA("acgtacctacgt"),
      makeDNA("acgtacgtacgt"),
      makeDNA("acgtacttacgt"),
    };

    final byte[][] expectedSequencesWithoutN = {
      makeDNA("acgtacgtacgt"),
      makeDNA("acgtacgtacgt"),
      makeDNA("acgtacgtacgt"),
      makeDNA("acgtacgtacgt"),
    };

    final String[] expectedSequenceNamesLeftWithN = {
      "test-exp-left-N->A", "test-exp-left-N->C", "test-exp-left-N->G", "test-exp-left-N->T",
    };

    final String[] expectedSequenceNamesRightWithN = {
      "test-exp-right-N->A", "test-exp-right-N->C", "test-exp-right-N->G", "test-exp-right-N->T",
    };

    final File dir = FileUtils.createTempDir("sdffilter", "testpair");
    try {
      final File leftInDir = new File(dir, "left-in");
      final File rightInDir = new File(dir, "right-in");
      final File outDir = new File(dir, "out");
      assertTrue(leftInDir.mkdir());
      assertTrue(rightInDir.mkdir());

      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(testInputStringWithN));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, leftInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();

      al = new ArrayList<>();
      al.add(createStream(testInputStringWithoutN));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, rightInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[] {"-a", "-l", leftInDir.getPath(), "-r", rightInDir.getPath(), "-o", outDir.getPath(), "-n", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());

      checkExpected(new File(outDir, "left"), expectedSequencesWithN, expectedSequenceNamesLeftWithN);
      checkExpected(new File(outDir, "right"), expectedSequencesWithoutN, expectedSequenceNamesLeftWithN);

      assertTrue(FileHelper.deleteAll(outDir));

      SdfFilter.mainInit(new String[] {"-a", "-l", leftInDir.getPath(), "-r", rightInDir.getPath(), "-o", outDir.getPath()}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());

      SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "left"));
      try {
        assertEquals("test", dsr.name(0));
        assertEquals(Arrays.toString(makeDNA(withN)), Arrays.toString(dsr.read(0)));
      } finally {
        dsr.close();
      }
      dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "right"));
      try {
        assertEquals("test", dsr.name(0));
        assertEquals(Arrays.toString(makeDNA(withoutN)), Arrays.toString(dsr.read(0)));
      } finally {
        dsr.close();
      }

      assertTrue(FileHelper.deleteAll(leftInDir));
      assertTrue(FileHelper.deleteAll(rightInDir));
      assertTrue(FileHelper.deleteAll(outDir));

      al = new ArrayList<>();
      al.add(createStream(testInputStringWithoutN));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, leftInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();

      al = new ArrayList<>();
      al.add(createStream(testInputStringWithN));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, rightInDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      SdfFilter.mainInit(new String[] {"-a", "-l", leftInDir.getPath(), "-r", rightInDir.getPath(), "-o", outDir.getPath(), "-n", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());

      checkExpected(new File(outDir, "left"), expectedSequencesWithoutN, expectedSequenceNamesRightWithN);
      checkExpected(new File(outDir, "right"), expectedSequencesWithN, expectedSequenceNamesRightWithN);
    } finally {
      assertTrue("Unable to delete " + dir + ": ", FileHelper.deleteAll(dir));
    }
  }

  public void testHelp() {
    TestCFlags.check(SdfFilter.getCFlags(System.out, System.err),
        "rtg sdffilter [OPTION]... -o DIR",
        "number of Ns to expand",
        "output base directory (must be empty if present)",
        "print help on command-line flag usage",
        "input SDF directory",
        "left input SDF directory",
        "right input SDF directory",
        "reverse complement all sequences",
        "alternate expansion for CG data"
    );
  }
}
