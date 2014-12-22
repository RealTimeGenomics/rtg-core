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
package com.rtg.index.hash;

import static com.rtg.mode.SequenceMode.BIDIRECTIONAL;
import static com.rtg.mode.SequenceMode.PROTEIN;
import static com.rtg.mode.SequenceMode.TRANSLATED;
import static com.rtg.mode.SequenceMode.UNIDIRECTIONAL;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractHashLoopTest extends TestCase {

  private File mDir;
  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() throws Exception {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  protected SequencesReader getReaderDNA(final String str) throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(new ByteArrayInputStream(str.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
        new DNAFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(mDir);
  }

  private SequencesReader getReaderProtein(final String str) throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(new ByteArrayInputStream(str.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al,
        new ProteinFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, mDir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();
    return SequencesReaderFactory.createDefaultSequencesReader(mDir);
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1() throws IOException {
    //tests stop codon TGA
    final String str0 = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG";
    assertEquals(64, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 1, 1, UNIDIRECTIONAL, 64, 1);
      checkCount(sr, 1, 1, BIDIRECTIONAL, 128, 2);
      checkCount(sr, 1, 1, TRANSLATED, 109, 6);

      checkCount(sr, 2, 1, UNIDIRECTIONAL, 63, 1);
      checkCount(sr, 2, 1, BIDIRECTIONAL, 126, 2);
      checkCount(sr, 2, 1, TRANSLATED, 89, 6);

      checkCount(sr, 2, 2, UNIDIRECTIONAL, 32, 1);
      checkCount(sr, 2, 2, BIDIRECTIONAL, 64, 2);
      checkCount(sr, 2, 2, TRANSLATED, 45, 6);

      checkCount(sr, 5, 1, UNIDIRECTIONAL, 60, 1);
      checkCount(sr, 5, 1, BIDIRECTIONAL, 120, 2);
      checkCount(sr, 5, 1, TRANSLATED, 50, 6);

      checkCount(sr, 5, 5, UNIDIRECTIONAL, 12, 1);
      checkCount(sr, 5, 5, BIDIRECTIONAL, 24, 2);
      checkCount(sr, 5, 5, TRANSLATED, 12, 6);

      checkCount(sr, 12, 12, TRANSLATED, 3, 6);

      checkCount(sr, 32, 32, UNIDIRECTIONAL, 2, 1);
      checkCount(sr, 32, 32, BIDIRECTIONAL, 4, 2);
    }
  }

  public final void testRelativelyPrimeLengthAndStep() throws IOException {
    final String str0 = "TAACCCCTCTC";
    assertEquals(11, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 11, 2, UNIDIRECTIONAL, 1, 1);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1n1() throws IOException {
    //check unknown codons N
    final String str0 = "CGTTGACGT";
    assertEquals(9, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 3, 1, TRANSLATED, 1, 6);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1n2() throws IOException {
    //check unknown codons N
    final String str0 = "CGTTAACGT";
    assertEquals(9, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 3, 1, TRANSLATED, 1, 6);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1n3() throws IOException {
    //check unknown codons N
    final String str0 = "CGTTAGCGT";
    assertEquals(9, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 3, 1, TRANSLATED, 1, 6);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1nx() throws IOException {
    //check unknown codons N
    final String str0 = "CGTTACCGT";
    assertEquals(9, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 3, 1, TRANSLATED, 2, 6);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test1a() throws IOException {
    final String str0 = "ARNDCQEGARNDCQEGARNDCQEGARNDCQEGARNDCQEGARNDCQEGARNDCQEGARNDCQEG";
    assertEquals(64, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderProtein(str)) {
      checkCount(sr, 2, 1, PROTEIN, 63, 1);
      checkCount(sr, 1, 1, PROTEIN, 64, 1);
      checkCount(sr, 2, 2, PROTEIN, 32, 1);
      checkCount(sr, 5, 1, PROTEIN, 60, 1);
      checkCount(sr, 5, 5, PROTEIN, 12, 1);
      checkCount(sr, 12, 12, PROTEIN, 5, 1);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test2() throws IOException {
    //tests stop codon TGA
    final String str0 = "ACTGA";
    assertEquals(5, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount(sr, 1, 1, UNIDIRECTIONAL, 5, 1);
      checkCount(sr, 1, 1, BIDIRECTIONAL, 10, 2);
      checkCount(sr, 1, 1, TRANSLATED, 5, 6);

      checkCount(sr, 2, 1, UNIDIRECTIONAL, 4, 1);
      checkCount(sr, 2, 1, BIDIRECTIONAL, 8, 2);
      checkCount(sr, 2, 1, TRANSLATED, 0, 6);

      checkCount(sr, 2, 2, UNIDIRECTIONAL, 2, 1);
      checkCount(sr, 2, 2, BIDIRECTIONAL, 4, 2);
      checkCount(sr, 2, 2, TRANSLATED, 0, 6);

      checkCount(sr, 5, 1, UNIDIRECTIONAL, 1, 1);
      checkCount(sr, 5, 1, BIDIRECTIONAL, 2, 2);
      checkCount(sr, 5, 1, TRANSLATED, 0, 6);

      checkCount(sr, 5, 5, UNIDIRECTIONAL, 1, 1);
      checkCount(sr, 5, 5, BIDIRECTIONAL, 2, 2);
      checkCount(sr, 5, 5, TRANSLATED, 0, 6);

      checkCount(sr, 12, 12, TRANSLATED, 0, 6);

      checkCount(sr, 32, 32, UNIDIRECTIONAL, 0, 1);
      checkCount(sr, 32, 32, BIDIRECTIONAL, 0, 2);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test2a() throws IOException {
    final String str0 = "ARNDC";
    assertEquals(5, str0.length());
    final String str = ""
      + ">x" + StringUtils.LS
      + str0 + StringUtils.LS;
    try (SequencesReader sr = getReaderProtein(str)) {
      checkCount(sr, 1, 1, PROTEIN, 5, 1);
      checkCount(sr, 2, 1, PROTEIN, 4, 1);
      checkCount(sr, 2, 2, PROTEIN, 2, 1);
      checkCount(sr, 5, 1, PROTEIN, 1, 1);
      checkCount(sr, 5, 5, PROTEIN, 1, 1);
      checkCount(sr, 12, 12, PROTEIN, 0, 1);
    }
  }

  protected abstract void checkCount(final SequencesReader sr, final int windowSize, final int stepSize, final SequenceMode mode, final int expected, final int maxId) throws IOException;

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test3() throws IOException {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ACT" + StringUtils.LS
      + ">x2" + StringUtils.LS
      + "ACTG" + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount2(sr, 1, 1, TRANSLATED, 6, 12);
      checkCount2(sr, 1, 1, UNIDIRECTIONAL, 7, 2);
      checkCount2(sr, 1, 1, BIDIRECTIONAL, 14, 4);

      checkCount2(sr, 2, 1, UNIDIRECTIONAL, 5, 2);
      checkCount2(sr, 2, 1, BIDIRECTIONAL, 10, 4);
      checkCount2(sr, 2, 1, TRANSLATED, 0, 12);

      checkCount2(sr, 2, 2, UNIDIRECTIONAL, 3, 2);
      checkCount2(sr, 2, 2, BIDIRECTIONAL, 6, 4);
      checkCount2(sr, 2, 2, TRANSLATED, 0, 12);

      checkCount2(sr, 5, 1, UNIDIRECTIONAL, 0, 2);
      checkCount2(sr, 5, 1, BIDIRECTIONAL, 0, 4);
      checkCount2(sr, 5, 1, TRANSLATED, 0, 12);
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test3a() throws IOException {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ARN" + StringUtils.LS
      + ">x2" + StringUtils.LS
      + "DCQE" + StringUtils.LS;
    try (SequencesReader sr = getReaderProtein(str)) {
      checkCount2(sr, 1, 1, PROTEIN, 7, 2);
      checkCount2(sr, 2, 1, PROTEIN, 5, 2);
      checkCount2(sr, 2, 2, PROTEIN, 3, 2);
      checkCount2(sr, 5, 1, PROTEIN, 0, 2);
    }
  }

  protected abstract void checkCount2(final SequencesReader sr, final int windowSize, final int stepSize, final SequenceMode mode, final int expected, final int maxId) throws IOException;

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public void test4() throws IOException {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ACGTACG" + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount3(sr, 1, 1, UNIDIRECTIONAL,
        new long[]{0, 1, 2, 3, 0, 1, 2},
        new int[]{0, 0, 0, 0, 0, 0, 0}
      );
      checkCount3(sr, 1, 1, BIDIRECTIONAL,
        new long[]{0, 1, 2, 3, 0, 1, 2,
          1, 2, 3, 0, 1, 2, 3},
        new int[]{0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1}
      );
      checkCount3(sr, 1, 1, TRANSLATED,
        new long[]{16, 18, 1, 16, 19, 1, 16, 19, 1, 18},
        new int[]{0, 0, 1, 1, 2, 3, 3, 4, 4, 5}
      );
      checkCount3(sr, 2, 1, UNIDIRECTIONAL,
        new long[]{1, 6, 11, 12, 1, 6},
        new int[]{0, 0, 0, 0, 0, 0}
      );
      checkCount3(sr, 2, 1, BIDIRECTIONAL,
        new long[]{1, 6, 11, 12, 1, 6,
          6, 11, 12, 1, 6, 11},
        new int[]{0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1}
      );
      checkCount3(sr, 2, 1, TRANSLATED,
        new long[]{530, 48, 48, 609},
        new int[]{0, 1, 3, 4}
      );
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void test4a() throws IOException {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ARNDCQE" + StringUtils.LS;
    try (SequencesReader sr = getReaderProtein(str)) {
      checkCount3(sr, 1, 1, PROTEIN,
        new long[]{0, 1, 2, 3, 4, 5, 6},
        new int[]{0, 0, 0, 0, 0, 0, 0}
      );
      checkCount3(sr, 2, 1, PROTEIN,
        new long[]{1, 34, 67, 100, 133, 166},
        new int[]{0, 0, 0, 0, 0, 0}
      );
    }
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public void test5() throws IOException {
    Diagnostic.setLogStream();
    final String str = ""
      + ">x" + StringUtils.LS
      + "ACGT" + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount3(sr, 2, 2, SequenceMode.BIDIRECTIONAL,
        new long[]{1, 11, 1, 11},
        new int[]{0, 0, 1, 1}
      );
    }
  }

  protected abstract void checkCount3(
      final SequencesReader sr, final int windowSize, final int stepSize,
      final SequenceMode mode, final long[] expectedL, final int[] expectedI)
  throws IOException;

  protected abstract void getLongLoop() throws IOException;

  public final void testLongError() throws IOException {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    try (PrintStream pr = new PrintStream(ba)) {
      Diagnostic.setLogStream(pr);
      try {
        getLongLoop();
      } catch (final SlimException e) {
        //expected
      } finally {
        Diagnostic.setLogStream();
      }
    } finally {
      ba.close();
    }
    final String str = ba.toString();
    assertTrue(str.contains("There is a sequence which is too long to process. Its length is \"2147483647\" bytes. See the SDF output for the name of the sequence."));
    //System.err.println(str);
  }
}

