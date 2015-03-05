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

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class SamAssistanceCgTest extends TestCase {

  private static final Locale DEFAULT = Locale.getDefault();

  static SAMRecord makeSamRecordCG(final String readString, final String cigar, final String superCigar, final String readDelta) {
    final SAMRecord sam = SamAssistanceCgLegacyTest.makeSamRecord(readString, cigar);
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, superCigar);
    if (readDelta.length() > 0) {
      sam.setAttribute(SamUtils.CG_READ_DELTA, readDelta);
    }
    return sam;
  }

  static final String TEMPLATE1 = SamAssistanceCgLegacyTest.TEMPLATE1;

  public void testCGNoOverlap() throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCg();
    final SAMRecord sam = makeSamRecordCG("acgtacgt".toUpperCase(DEFAULT), "7X", "7X", "ACGTACG");
    final String[] actual = sama.samToReads(sam, TEMPLATE1, TEMPLATE1.getBytes(), 1, false);
    assertEquals(1, actual.length);
    assertEquals(" ACGTACG", actual[0]);
  }

  void check2Reads(final String line1, final String line2, final String read, final String cigar, final String superCigar, final String readDelta, boolean displayDots) throws BadSuperCigarException {
    check2Reads(line1, line2, read, cigar, superCigar, readDelta, displayDots, TEMPLATE1);
  }
  // for testing the complex CG cases, which are displayed on two lines.
  void check2Reads(final String line1, final String line2, final String read, final String cigar, final String superCigar, final String readDelta, boolean displayDots, String template) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCg();
    final SAMRecord sam = makeSamRecordCG(read.replaceAll(" ", "").toUpperCase(DEFAULT), cigar, superCigar, readDelta);
    final String[] actual = sama.samToReads(sam, template, DnaUtils.encodeString(template), 0, displayDots);
    //System.out.println(Arrays.toString(actual));
    assertEquals(2, actual.length);
    assertEquals(line1, actual[0]);
    assertEquals(line2, actual[1]);
  }

  public void testCG2AllEqual() throws BadSuperCigarException {
    check2Reads(
        "..........",
        "        ....",
        TEMPLATE1, "12=", "10=2B4=", "", true
    );

    check2Reads(
        "ACGTAACCGG",
        "        GGTT",
        TEMPLATE1, "12=", "10=2B4=", "", false
    );
  }

  public void testCG2SnpInBothOverlap() throws BadSuperCigarException {
    // a SNP in both parts of the overlap region
    check2Reads(
        ".......T.A",
        "        C..A",
        "acgtaacTgAtA", "7=1X1=1X1=1X", "7=1X1=1X2B1X2=1X", "TACA", true
    );
  }

  public void testCG2SnpInSecondOverlap() throws BadSuperCigarException {
    // a SNP in the second part of the overlap region
    check2Reads(
        ".......T..",
        "        C..A",
        "acgtaacTggtA", "7=1X3=1X", "7=1X2=2B1X2=1X", "TCA", true
    );
  }

  public void testCG2Deletion() throws BadSuperCigarException {
    // a deletion in the read
    check2Reads(
        ".....--AAA",
        "        .C..",
        "acgtaAAAtt", "5=2D3X2=", "5=2D3X2B1=1X2=", "AAAC", true
    );
  }

  public void testCG2NsAndLotsOfSnps() throws BadSuperCigarException {
    // an N gap in the read, an overlap of 3, and SNPs all through the both overlaps!
    check2Reads(
        "...   ..ACT",
        "        TCAG",
        "acg   ccACTG", "3=3N2=4X", "3=3N2=3X3B4X", "ACTTCAG", true
    );
  }

  public void testCG2InitialInserts() throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCg();
    final SAMRecord sam = makeSamRecordCG(TEMPLATE1.replaceAll(" ", "").toUpperCase(DEFAULT), "12=", "10=2B4=", "");
    final String template = "a__cgt_aaccggtt";
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true);
    assertEquals(2, actual.length);
    assertEquals(".__..._......", actual[0]);
    assertEquals("           ....", actual[1]);
  }

  private void checkOverlap2(final String template, final String exp0, final String exp1) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCg();
    final SAMRecord sam = makeSamRecordCG("acgtaacTgAtA".toUpperCase(DEFAULT), "7=1X1=1X1=1X", "7=1X1=1X2B1X2=1X", "TACA");
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true);
    assertEquals(2, actual.length);
    assertEquals(exp0, actual[0]);
    assertEquals(exp1, actual[1]);
  }

  public void testCG2InsertsOnOverlap2() throws BadSuperCigarException {
    checkOverlap2("a__cgt_aacc_ggtt", ".__..._...T_.A", "            C..A");
    checkOverlap2("a__cgt_aacc___ggtt", ".__..._...T___.A", "              C..A");

    checkOverlap2("a__cgt_aaccg_gtt", ".__..._...T._A",
                                      "           C_..A");
    checkOverlap2("a__cgt_aaccg___gtt", ".__..._...T.___A", "           C___..A");

    checkOverlap2("a__cgt_aaccgg_tt", ".__..._...T.A", "           C._.A");
    checkOverlap2("a__cgt_aaccgg___tt", ".__..._...T.A", "           C.___.A");
  }

  private void checkOverlap1(final String template, final String exp0, final String exp1) throws BadSuperCigarException {
    final SamAssistance sama = new SamAssistanceCg();
    final SAMRecord sam = makeSamRecordCG("acgtaacTgAtA".toUpperCase(DEFAULT), "7=1X1=1X1=1X", "7=1X1=1B2X1=1X", "TCAA");
    final String[] actual = sama.samToReads(sam, template, template.getBytes(), 0, true);
    assertEquals(2, actual.length);
    assertEquals(exp0, actual[0]);
    assertEquals(exp1, actual[1]);
  }

  public void testCG2InsertsOnOverlap1() throws BadSuperCigarException {
    checkOverlap1("a__cgt_aacc_ggtt", ".__..._...T_.", "            CA.A");
    checkOverlap1("a__cgt_aacc___ggtt", ".__..._...T___.", "              CA.A");

    checkOverlap1("a__cgt_aaccg_gtt", ".__..._...T.", "           C_A.A");
    checkOverlap1("a__cgt_aaccg___gtt", ".__..._...T.", "           C___A.A");

    checkOverlap1("a__cgt_aaccgg_tt", ".__..._...T.", "           CA_.A");
    checkOverlap1("a__cgt_aaccgg___tt", ".__..._...T.", "           CA___.A");
  }

  public void testCGInserts() throws BadSuperCigarException {
    check2Reads("....."
              , "   ....................      ......A.A."
              , "AAGGATCACAGAGGTTCCAGAGTGAACTGAGAG", "23=6N6=1I1=1X1=", "5=2B20=6N6=1I1=1X1=", "AA", true
              , "aaggatcacagaggttccagagtcagggagaactg_gc");
  }

}
