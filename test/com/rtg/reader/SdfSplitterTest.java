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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.mode.DNA;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.Residue;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test for PrereadSplitter
 *
 */
public class SdfSplitterTest extends AbstractCliTest {



  public static Test suite() {
    return new TestSuite(SdfSplitterTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.runAndWait(suite());
  }

  @Override
  protected AbstractCli getCli() {
    return new SdfSplitter();
  }

  private InputStream createStream(final String data) {
    return new ByteArrayInputStream(data.getBytes());
  }

  private static final String TEST_INPUT_STRING = ">test\n"
      + "acgtgtgtgtcttagggctcactggtcatgca\n"
      + ">bob the buuilder\n"
      + "tagttcagcatcgatca\n"
      + ">hobos r us\n"
      + "accccaccccacaaacccaa";

  private static final String TEST_INPUT_STRING_2 = ">test2\n"
      + "acgtgtgtgtcttagggctcactggtcatgca\n"
      + ">bob2 the buuilder\n"
      + "tagttcagcatcgatca\n"
      + ">hobos3 r us2\n"
      + "accccaccccacaaacccaa"
      + ">hobos4 r us3\n"
      + "accccaccccacaaacccaa";

  private static final Residue[][] EXPECTED_SEQUENCES;
  private static final Residue[][] EXPECTED_SEQUENCES_2;
  static {
    Residue a = DNA.A;
    Residue c = DNA.C;
    Residue g = DNA.G;
    Residue t = DNA.T;
    EXPECTED_SEQUENCES = new Residue[][] {
      //test
      {a, c, g, t, g, t,
      g, t, g, t, g, t,
      g, t, c, t, t, a,
      t, a, g, g, g, c,
      g, c, t, c, a, c,
      a, c, t, g, g, t,
      g, t, c, a, t, g,
      t, g, c, a},
      //bob the builder
      {t, a, g, t, t, c,
      t, c, a, g, c, a,
      c, a, t, c, g, a,
      g, a, t, c, a},
      //hobos r us
      {a, c, c, c, c, a,
      c, a, c, c, c, c,
      c, c, a, c, a, a,
      a, a, a, c, c, c,
      c, c, a, a}
    };

    EXPECTED_SEQUENCES_2 = new Residue[][] {
      //test2
      {a, c, g, t, g, t,
      g, t, g, t, g, t,
      g, t, c, t, t, a,
      t, a, g, g, g, c,
      g, c, t, c, a, c,
      a, c, t, g, g, t,
      g, t, c, a, t, g,
      t, g, c, a},
      //bob2 the builder
      {t, a, g, t, t, c,
      t, c, a, g, c, a,
      c, a, t, c, g, a,
      g, a, t, c, a},
      //hobos.1 r us
      {a, c, c, c, c, a,
      c, a, c, c, c, c,
      c, c, a, c, a, a,
      a, a, a, c, c, c,
      c, c, a, a},
      //hobos.2 r us
      {a, c, c, c, c, a,
      c, a, c, c, c, c,
      c, c, a, c, a, a,
      a, a, a, c, c, c,
      c, c, a, a}
    };
  }

  private static final String[] EXPECTED_SEQUENCE_NAMES = {
    "test", "bob", "hobos"
  };
  private static final String[] EXPECTED_SEQUENCE_NAMES_2 = {
    "test2", "bob2", "hobos3", "hobos4"
  };


  public void testNotSplitting() throws IOException {
    Diagnostic.setLogStream();
    try (TestDirectory inDir = new TestDirectory("prereadsplitter-test.in");
         TestDirectory outDir = new TestDirectory("prereadsplitter-testout")) {
      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
              new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      ArrayList<File> inDirs = new ArrayList<>();
      inDirs.add(inDir);
      SdfSplitter.split(inDirs, outDir, 100, true, false, false);
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "000000"))) {
        assertEquals(EXPECTED_SEQUENCES.length, dsr.numberSequences());
        for (int j = 0; j < EXPECTED_SEQUENCES.length; j++) {
          assertEquals("j: " + j, EXPECTED_SEQUENCE_NAMES[j], dsr.name(j));
        }
      }
    }
  }

  public void testSplitting() throws IOException {
    Diagnostic.setLogStream();
    try (TestDirectory inDir = new TestDirectory("prereadsplitter-test.in");
         TestDirectory outDir = new TestDirectory("prereadsplitter-testout")) {
      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
              new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      ArrayList<File> inDirs = new ArrayList<>();
      inDirs.add(inDir);
      SdfSplitter.split(inDirs, outDir, 1, true, false, false);
      for (int j = 0; j < EXPECTED_SEQUENCES.length; j++) {
        SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "00000" + j));
        SequencesIterator it = dsr.iterator();
        assertEquals("j: " + j, 1, dsr.numberSequences());
        try {
          assertTrue(it.nextSequence());
          assertEquals("j: " + j, EXPECTED_SEQUENCE_NAMES[j], it.currentName());
        } finally {
          dsr.close();
        }
      }
    }
  }

  public void testMerging() throws IOException {
    try (TestDirectory inDir = new TestDirectory("prereadsplitter-test.in");
         TestDirectory inDir2 = new TestDirectory("prereadsplitter-test.in2");
         TestDirectory outDir = new TestDirectory("prereadsplitter-testout")) {
      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      ArrayList<File> inDirs = new ArrayList<>();
      inDirs.add(inDir);

      al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING_2));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, inDir2, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      inDirs.add(inDir2);

      SdfSplitter.split(inDirs, outDir, 100, false, false, false);
      try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(new File(outDir, "000000"))) {
        SequencesIterator it = dsr.iterator();
        for (int j = 0; j < EXPECTED_SEQUENCES.length; j++) {
          assertTrue("j: " + j, it.nextSequence());
          assertEquals("j: " + j, EXPECTED_SEQUENCE_NAMES[j], it.currentName());
        }
        for (int j = 0; j < EXPECTED_SEQUENCES_2.length; j++) {
          assertTrue("j2: " + j, it.nextSequence());
          assertEquals("j2: " + j, EXPECTED_SEQUENCE_NAMES_2[j], it.currentName());
        }
      }
    }
  }

  public void testDuplicates() throws IOException {
    Diagnostic.setLogStream();
    try (TestDirectory inDir = new TestDirectory("prereadsplitter-test.in");
         TestDirectory inDir2 = new TestDirectory("prereadsplitter-test.in2");
         TestDirectory tempDir = new TestDirectory("prereadsplitter-testout")) {
      final File outDir = new File(tempDir, "output");
      final File dupFile = new File(outDir.getAbsolutePath(), "duplicate-names.txt");
      ArrayList<InputStream> al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING));
      FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      SequencesWriter sw = new SequencesWriter(ds, inDir, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      al = new ArrayList<>();
      al.add(createStream(TEST_INPUT_STRING));
      ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      sw = new SequencesWriter(ds, inDir2, 20, PrereadType.UNKNOWN, false);
      sw.processSequences();
      final String err = checkMainInitWarn("-n", "100", "-o", outDir.getAbsolutePath(), inDir.getAbsolutePath(), inDir2.getAbsolutePath());
      assertTrue(outDir.exists());
      assertTrue(dupFile.exists());
      assertTrue(err.contains("Duplicate Sequence Names in Input"));
      final String[] expected = {"bob" + StringUtils.LS, "test" + StringUtils.LS, "hobos" + StringUtils.LS};
      TestUtils.containsAll(FileUtils.fileToString(dupFile), expected);
    }
  }

  public void testValidArgs() throws IOException {
    try (TestDirectory inDir = new TestDirectory("prereadsplitter-test.in");
         TestDirectory outDir = new TestDirectory("prereadsplitter-testout")) {
      assertTrue(outDir.delete());
      checkHandleFlagsOut("-n", "1000000", "-o", outDir.getPath(), inDir.getPath());
    }
  }

  public void testHelp() {
    checkHelp("sdfsplit [OPTION]... -n INT -o DIR",
              "[OPTION]... -n INT -o DIR -I FILE",
              "number of reads per output",
              "output base directory",
              "input SDF",
              "file containing a list of SDFs (1 per line)");
  }
}
