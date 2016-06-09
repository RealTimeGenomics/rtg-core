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
package com.rtg.visualization;

import java.util.Locale;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class SamAssistanceCgLegacyTest extends TestCase {

  private static final Locale DEFAULT = Locale.getDefault();

  static SAMRecord makeSamRecord(final String readString, final String cigar) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString(readString);
    sam.setCigarString(cigar);
    sam.setBaseQualities(new byte[readString.length()]);
    return sam;
  }

  public void testSplitCigar() {
    checkSplitCigar("10I5=", 10, "10I", "5=");
    checkSplitCigar("9X5=2X", 10, "9X1=", "4=2X");
    try {
      checkSplitCigar("5=3M", 10, "5=3M", null);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("bad CG Cigar=5=3M split=10", e.getMessage());
    }
    checkSplitCigar("8=3X", 10, "8=2X", "1X");
    checkSplitCigar("8=4N3X", 10, "8=4N2X", "1X");
    checkSplitCigar("8=4D3X", 10, "8=4D2X", "1X");
    checkSplitCigar("8S3X", 10, "8S2X", "1X");
    checkSplitCigar("8=3X5I2=", 10, "8=2X", "1X5I2=");
    checkSplitCigar("1X5I2=", 5, "1X4I", "1I2=");

    try {
      checkSplitCigar("10I5=", 0, null, "10I5=");
    } catch (final IllegalArgumentException e) {
      assertEquals("bad CG split=0", e.getMessage());
    }

    try {
      checkSplitCigar("8Z4D3X", 10, "8=4D2X", "1X");
    } catch (final IllegalStateException e) {
      assertEquals("Invalid cigar : 8Z4D3X", e.getMessage());
    }
  }

  private void checkSplitCigar(final String cigar, final int split, final String exp0, final String exp1) {
    final String[] actual = SamAssistanceCgLegacy.splitCigar(cigar, split);
    assertEquals(2, actual.length);
    assertEquals(exp0, actual[0]);
    assertEquals(exp1, actual[1]);
  }

  void check1Read(final String expected, final String read, final String cigar, final int[] inserts) throws BadSuperCigarException {
    check1Read(expected, read, cigar, inserts, 0, true);
  }

  void check1Read(final String expected, final String read, final String cigar, final int[] inserts, final int readStart, final boolean displayDots) throws BadSuperCigarException {
    check1Read(expected, read, cigar, SamAssistanceSimpleTest.insertString(inserts), readStart, displayDots);
  }

  private void check1Read(String expected, String read, String cigar, String template, int readStart, boolean displayDots) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecord(read, cigar);
    sam.setAlignmentStart(1);
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), readStart, displayDots, false);
    assertEquals(1, actual.length);
    assertEquals(expected, actual[0]);
  }

  public final void testCigarToReadsWithInserts() throws BadSuperCigarException {
    check1Read("ACGT", "ACGT", "4M", new int[] {0, 0, 0, 0});
    check1Read("A  C", "ACGT", "1M2N1M", new int[] {0, 0, 0, 0});

    check1Read("A--C", "ACGT", "1M2D1M", new int[] {0, 0, 0, 0});
    check1Read("....", "ACGT", "4=", "ACGT", 0, true);
    check1Read("ACGT", "ACGT", "4=", "NNNN", 0, true);
    check1Read("..N.", "ACNT", "4=", "ACGT", 0, true);
    check1Read("ACGT", "ACGT", "4X", new int[] {0, 0, 0, 0});

    //test "S"
    check1Read("CG.", "ACGT", "1S2X1=", "TTT", 0, true);
    check1Read("ACNT", "ACNT", "4X", new int[] {0, 0, 0, 0});

    check1Read("AC_GT", "ACGT", "4M", new int[] {0, 0, 1, 0});
    check1Read(".._..", "ACGT", "4=", "AC_GT", 0, true);
    check1Read("AC_GT", "ACGT", "4X", new int[] {0, 0, 1, 0});
    check1Read("ACG____T", "ACGT", "2M1I1M", new int[] {0, 0, 5, 0});

    check1Read("..GTGT", "ACGTGT", "2=1I1X1I1X", "AC_A_A", 0, true);

    check1Read(".-CG", "ACGT", "1=1D2X", "ACTT", 0, true);

    check1Read("   .-CGT", "ACGT", "1=1D3X", "___AAAAA", 3, true);

    check1Read("A-CG", "ACGT", "1=1D2X", new int[] {0, 0, 0, 0}, 0, false);
  }

  static SAMRecord makeSamRecordLegacyFormatCG(final String readString, final String cigar, final String gs, final String gc) {
    final SAMRecord sam = makeSamRecord(readString, cigar);
    sam.setAttribute("GS", gs);
    sam.setAttribute("GC", gc);
    sam.setAttribute("GQ", new byte[gs.length() / 2]);
    return sam;
  }

  static final String TEMPLATE1 = "acgtaaccggtt"; // length = 12

  public void testCGNoOverlapLegacy() throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecord("acgtacgt".toUpperCase(DEFAULT), "7X");
    final String[] actual = sama.samToReads(sam, TEMPLATE1, TEMPLATE1.getBytes(), 1, false, false);
    assertEquals(1, actual.length);
    assertEquals(" ACGTACG", actual[0]);
  }

  // for testing the complex CG cases, which are displayed on two lines.
  void check2ReadsLegacy(final String line1, final String line2, final String read, final String cigar, final String gs, final String gc) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecordLegacyFormatCG(read.replaceAll(" ", "").toUpperCase(DEFAULT), cigar, gs.toUpperCase(DEFAULT), gc.toUpperCase(DEFAULT));

    final String[] actual = sama.samToReads(sam, TEMPLATE1, TEMPLATE1.getBytes(), 0, true, false);
    assertEquals(2, actual.length);
    assertEquals(line1, actual[0]);
    assertEquals(line2, actual[1]);
  }

  public void testCG2AllEqualLegacy() throws BadSuperCigarException {
    check2ReadsLegacy(
        "..........",
        "        ....",
        TEMPLATE1, "12=", "gggg", "8S2G2S"
    );
  }

  public void testCG2SnpInBothOverlapLegacy() throws BadSuperCigarException {
    // a SNP in both parts of the overlap region
    check2ReadsLegacy(
        ".......T.A",
        "        C..A",
        "acgtaacTgAtA", "7=1X1=1X1=1X", "ggCg", "8S2G2S"
    );
  }

  public void testCG2SnpInSecondOverlapLegacy() throws BadSuperCigarException {
    // a SNP in the second part of the overlap region
    check2ReadsLegacy(
        ".......T..",
        "        C..A",
        "acgtaacTggtA", "7=1X3=1X", "ggCg", "8S2G2S"
    );
  }

  public void testCG2DeletionLegacy() throws BadSuperCigarException {
    // a deletion in the read
    check2ReadsLegacy(
        ".....--AAA",
        "        .C..",
        "acgtaAAAtt", "5=2D3X2=", "gggC", "6S2G2S"
    );
  }

  public void testCG2NsAndLotsOfSnpsLegacy() throws BadSuperCigarException {
    // an N gap in the read, an overlap of 3, and SNPs all through the both overlaps!
    check2ReadsLegacy(
        "...   ..AC.",
        "        TCAG",
        "acg   ccACTG", "3=3N2=4X", "xxxTCA", "5S3G1S"
    );
  }

  public void testCG2InitialInsertsLegacy() throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecordLegacyFormatCG(TEMPLATE1.replaceAll(" ", "").toUpperCase(DEFAULT), "12=", "gggg".toUpperCase(DEFAULT), "8S2G2S".toUpperCase(DEFAULT));
    final String template = "a__cgt_aaccggtt";
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true, false);
    assertEquals(2, actual.length);
    assertEquals(".__..._......", actual[0]);
    assertEquals("           ....", actual[1]);
  }

  private void checkOverlap2Legacy(final String template, final String exp0, final String exp1) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecordLegacyFormatCG("acgtaacTgAtA".toUpperCase(DEFAULT), "7=1X1=1X1=1X", "ggCg".toUpperCase(DEFAULT), "8S2G2S".toUpperCase(DEFAULT));
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true, false);
    assertEquals(2, actual.length);
    assertEquals(exp0, actual[0]);
    assertEquals(exp1, actual[1]);
  }

  public void testCG2InsertsOnOverlap2Legacy() throws BadSuperCigarException {
    checkOverlap2Legacy("a__cgt_aacc_ggtt", ".__..._...T_.A", "            C..A");
    checkOverlap2Legacy("a__cgt_aacc___ggtt", ".__..._...T___.A", "              C..A");

    checkOverlap2Legacy("a__cgt_aaccg_gtt", ".__..._...T._A", "           C_..A");
    checkOverlap2Legacy("a__cgt_aaccg___gtt", ".__..._...T.___A", "           C___..A");

    checkOverlap2Legacy("a__cgt_aaccgg_tt", ".__..._...T.A", "           C._.A");
    checkOverlap2Legacy("a__cgt_aaccgg___tt", ".__..._...T.A", "           C.___.A");
  }

  private void checkOverlap1Legacy(final String template, final String exp0, final String exp1) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCgLegacy();
    final SAMRecord sam = makeSamRecordLegacyFormatCG("acgtaacTgAtA".toUpperCase(DEFAULT), "7=1X1=1X1=1X", "gC".toUpperCase(DEFAULT), "8S1G3S".toUpperCase(DEFAULT));
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true, false);
    assertEquals(2, actual.length);
    assertEquals(exp0, actual[0]);
    assertEquals(exp1, actual[1]);
  }

  public void testCG2InsertsOnOverlap1Legacy() throws BadSuperCigarException {
    checkOverlap1Legacy("a__cgt_aacc_ggtt", ".__..._...T_.", "            CA.A");
    checkOverlap1Legacy("a__cgt_aacc___ggtt", ".__..._...T___.", "              CA.A");

    checkOverlap1Legacy("a__cgt_aaccg_gtt", ".__..._...T.", "           C_A.A");
    checkOverlap1Legacy("a__cgt_aaccg___gtt", ".__..._...T.", "           C___A.A");

    checkOverlap1Legacy("a__cgt_aaccgg_tt", ".__..._...T.", "           CA_.A");
    checkOverlap1Legacy("a__cgt_aaccgg___tt", ".__..._...T.", "           CA___.A");
  }
}
