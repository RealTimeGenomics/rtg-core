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
package com.rtg.alignment;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;


/**
 * Tests corresponding class
 */
public class EditDistanceFactoryTest extends TestCase {

  public void testCreate() throws Exception {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    Diagnostic.setLogStream();
    try {
      EditDistanceFactory.createEditDistance(params, null, null);
      fail();
    } catch (final IllegalArgumentException iae) {
      assertEquals("reader1 cannot be null", iae.getMessage());
    }

    final File tmpDir = FileUtils.createTempDir("edFact", "test");
    try {
      final SequencesReader shortreader = ReaderTestUtils.getReaderDNA(">read1" + LS + "acgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgact" + LS, new File(tmpDir, "short"), null);

      final File badCGReadsDir = new File(tmpDir, "badcg");
      final SequencesReader cgreader = ReaderTestUtils.getReaderDNAFastqCG("", badCGReadsDir, PrereadArm.LEFT);
      try {
        EditDistanceFactory.createEditDistance(params, cgreader, null);
        fail();
      } catch (final IllegalArgumentException iae) {
        assertEquals("CG requires paired data", iae.getMessage());
      }
      try {
        EditDistanceFactory.createEditDistance(params, cgreader, shortreader);
        fail();
      } catch (final IllegalArgumentException iae) {
        // Expected
      }
      try {
        EditDistanceFactory.createEditDistance(params, PrereadType.CG, 10, 20);
        fail();
      } catch (final IllegalArgumentException iae) {
        assertEquals("CG data only supports " + CgGotohEditDistance.CG_RAW_READ_LENGTH + " length reads", iae.getMessage());
      }

      final File leftCGReadsDir = new File(tmpDir, "leftcg");
      final SequencesReader cgreader2 = ReaderTestUtils.getReaderDNAFastqCG("@s1" + LS + "tagacaaatgggattacaagccacaggagtgaata" + "\n+\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n", leftCGReadsDir, PrereadArm.LEFT);

      EditDistance e = EditDistanceFactory.createEditDistance(params, cgreader2, cgreader2);
      assertTrue(e instanceof RcEditDistance);
      /*
      assertTrue(e instanceof PrioritisedEditDistance);
      final PrioritisedEditDistance ee = (PrioritisedEditDistance) e;
      assertTrue(ee.mEds[0] instanceof CgPatternAligner);
      assertTrue(ee.mEds[1] instanceof LoopingEditDistance);
      assertTrue(((LoopingEditDistance) ee.mEds[1]).mEd instanceof CgEditDistance);
      */

      e = EditDistanceFactory.createEditDistance(params, shortreader, null);
      assertTrue(e instanceof SoftClipperOmni);

      //final SequencesReader longreader = ReaderTestUtils.getReaderDNA(">read1" + LS + "acgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgacgtgactg" + LS, new File(tmpDir, "long"));
      //e = EditDistanceFactory.createEditDistance(shortreader, longreader, 1);
      //assertTrue(e instanceof LongReadEditDistance);
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testCalculateLowerBoundValue() {
    Diagnostic.setLogStream();
    final MaxShiftFactor msf = new MaxShiftFactor(0.1);
    assertEquals(3, EditDistanceFactory.calculateLowerBoundValue(36, msf));
    assertEquals(4, EditDistanceFactory.calculateLowerBoundValue(100, msf));
    assertEquals(3, EditDistanceFactory.calculateLowerBoundValue(1000, msf));
    assertEquals(3, EditDistanceFactory.calculateLowerBoundValue(2000, msf));

  }

  public void testExample1() {
    // This one failed due to approximating modulo with & not working for non-powers-of-two. Now fixed.
    checkEarlyTerminationDoesntLie("TAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAACGTGTCATCTAGCATTAGGTATGTCTCCCGATGCTATCCCTCCCCCCTCCCCCGACCCCACCACAGTCCCCAGAGTGTGATATTCCCCTTCCTGTGTCCATGTCATCTCATTGTTCAATTCCTACCTATGAGTGAGAATATGCGGTGTTTGGTTTTTTGTTCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATAGTCCTTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTACTAGTTCTAGATCCGTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACAGTCCCACCAACAGTGTAAAAGTGTTCCTATTTCTCCACATCCTCTCCAGCACCTGTTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGATATCTCATAGTGATTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATTTCTTCATGTGTTTTTTGGCTGCATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATGTCCTTGGCCCACTTTTTGATGGGGTTGTTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGTAGGTTGCGAAAATTTTCTCCCATGTTGTAGGTTGCCTGTTCACTCTGATGGTAGTTTCTTTTGCTGTGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGGCTT",
                                   "TTTAGGGTACATGTGCTAATTGAGCAGGTTAGTTACATATGTATAACCATGTGCCATTGCTGGTGCGCTGCACCCACTAACGTGTCATCTAGCATTAGGTATGTCTCCCGATGCTATCCCTCCCCCCTCCCCCGACCCCACCACAGTCCCCAGAGTGTGATATTCCCTTTCCTGTGTCCATGTCATCTCATTGTTCAATTCCTACCTATGAGTGAGGAATATGCGGTGTTTGGTTTTTTGGGCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTAACAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCACATTTTCTTAATCCAGTCTATATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGTTATTGTGAATAGTGCCGCAATAAGCATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATAGTCCCTTGGGTATATACCCAGTAATGGGATGGCTGGGTGCAAATGGTATTACTAGTTCTAGATCCGTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACAGTCCCACCAACAGTGTAAAAGTGTGTTTATTTCTCCCCATCCTCTCCAGCACCTGTTGTTTCCTGACTTTTTAAAGATTGCCATTCTAACTGGTGTGAGATGCTATCTCTAAGATAGTGATTTTAATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATTTCTTCATTGTGTTTTTTGGCTGCATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATGTCCTTGGCCCACTTTTTGATGGGGTTGTTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGATATTAGCCCTTTGTCAGATGAGTAGGTTGCGAAAATTTTCTCCCATGTTGTAGGTTGCCTGTTCACTCCTGATGGTAGTTTCCTTTGCTGTGCAGAAGCTCTTTAGTTTAATTATCTCCCATT",
                                   0, 58, 18);
  }

  public void testFlailingTestExample2() {
    final String tmpl = "AGTAAGGGGGCTCCCTGCCCCCTGGGGGAACCTGACCCTAATTCTGACTCAGCAACTCACCATGAGCCATTATCCCGGGGATGCCAACAGCCCCCTTAGGACATGGCCTGAATGGACAGACACATAAGGCCTCGTGTGAGGAAGAAACACAGAAGCACCGGGCAGGGGGAAGCAGGCGCTCGAGGGGCTCACAGATGCACACACGACGAAGCCGTTCTTGTGGCTCAGCTCAAGGCTCAGCTGAACTGAAGAAGCAATTATAGCTTTTAATTCAAAAATAAGATTCCCAAGGATTGGTGCCATTTGCCCATTAACTTCACAGATTATACCTATTAGTGGTCAAGACAAACACGTGGCTCTCAGATCTAAACCTTCTTTAGGGACACTAAAGAAACTGAGGCTTGGGACAGACAGACACCAGGCAGAGCTGAGCAGGTAACTTAACAGCTGACACTGAGGGAGCAGAAAATGGCATTTCTCTCATTCTGCAGCAGACACGGGGCAGTGTGGGCCAGATGCCAAGTATGAGGCCAGGACACCAAGCCAGCCCTGGGCCCTCACTGCTCAGCCAGGACCCTTTTCTCCTGAACTAAACCATCTGGGTTTAGTTGAACTAAACCTGGGGCTGCCCGTGCAAAGGCATCTGCCACATTCAGTTCTCAGCAAGTTGGCAACACTGACAAGCATCCCTGCAAGAGAGGTGAAGACTCAAGTTCTATTTGTTTCTGGCATGAAAAGCGATGGTTTCCTTCTCATCACTCATAAAAAATAGCAGAAATGTGGCAATTGCGGGATCGATCCCACCGCACACCATACGTGTGCATATGCTACCATGCACATAGGCCTCAACCCAGGCATGGGAAATACCAGTAGGAGCGGCAAAGAGACCTTCCAATCTCTCTAGAAACTCCCTACAACTGAACCAGTGCTACTTCAGATCCCGCCGCACACCATACGCGTGCATATGCTACCATGCACATAGGACGTGAGGTCTGGCGTCAACTCAGGCATGGGAAATACCAGTAGGAGTGGCAAAGAGACCTTCCAATCTCTCTAGAAACTCCCTACAACTGGACCAGCGCTACTTCAGGTGACAGCCGCGTACATTCCTCAAATGATGACAACAGGCAGCCAGAATACCGC";
    final String read = "agtaaggggggctccctgccccctgggggaacctgaccctaattctgactcagcaactcaccatgagccattatcccggggatgccaacagcccccttaggacatggcctgaatggacagacacataaggcctcgtgtgaggaagaaacacagaagcaccgggcaggggggaagcaggcgctcgaggggctcacaagatgcacacacgacgaagccgttcttgtggctcagctcaaggctcagctgaactgaagaagcaattatagcttttaattcaaaagtaagattcccaaggattggtgccatttgcccattaacttcacagaatcgaacctattagtggtcaagacaaacacgtggctctcagatctaaaccttctttagggacactaaagaaactggaggcttgggacagacagacaccaggcagagctgagcaggtaacttaacagctggtactgagggagcagaaaatggcatttctctcattctgcagcagacacggggcagtgtgggccagatgccaagtatgaggccaggacaccaagccagccctgggccctcactgctcagccaggacccttttctcctgaactaaaccatctgggttttagttgacactaaacgtggggctgcccgtgcaaaggcatctgccacattcagttctcagcaagttggcaacactgacaagcatccctgcaagagaggtgaagactcaagttctatttgtttctggcatggaaaagcgatggttttttccttctcatcactcataaaaaatagcagaaatgtggcaattgcgggatcgatcccaccgcacaccatacgtgtgcatatgctaccatgcacataggcctcaacccaggcatgggaaataccagtaggagcggcaaagagaccttccaatctctctagaaacttcctacaactgaaccagtgctacttcagatcccgccgcacaccatacgcgtgcatatgctaccatgcacataggacgtga";
    // This one finds an alignment by shifting 142!!
    // Seeded aligner should be allowed to shift so far?  Lower bound estimator should allow some shifting
    checkEarlyTerminationDoesntLie(tmpl, read, 142, 100, 18);
//    final int maxShift = MaxShiftUtils.calculateDefaultMaxShift(read.length());
    checkEarlyTerminationDoesntLie(tmpl, read, 19, 100, 18);
    checkEarlyTerminationDoesntLie(tmpl, read, 18, 100, 18);
  }

  public void testFlailingTestExample3() {
    // This one finds a better alignment without shifting!!!
    checkEarlyTerminationDoesntLie("CTGGGACTACACCAGCTGCCACCATGCCTGGCTAATTTTTCGTATTTTTAGTAGAGACAGGGTTTCACTGAATTGGCCAGGATGGTCTTGATCTCCTCACCTTGTGATCCTCTTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCCTGAGCCAAGATACATATTTTTTAAAGGAAGAAAAATTTCAAAGGTACTCTGTTTGGTACAATAATCAAATATATAAATTGAGGAATAAAACATAACCACGAAACATATTTATAACTGCATATGGAAAATACAGAGGATAATTTTTTAAATAACATATTTTGAAAACCTTAACTAGGAATTTGAAAAGATCGCATTTGACAGGCCAGTATGAACACAACTTGAACGCAGCAAGACAGGTTCCCCATAAAAAAATCAAACACAGGGAAAATGAAACCACAAAGTTTCAATCTGCTCTGACCTTTGAAAAACTCAGCACAGACAGTGGCACTTAGGACCAACAGCAGGAGATCCCTAATCCCATCACCATGGCGATAGGGCATAAACATTCCAGGGTGAAGTCACAATCCACATTGGGAGGTCCAACTGCTGCCAGGCAGACAGGTGTGCCTTTACATGTACAGGAAGGTCATTGAAGGCTCAGTGTTTTGTTTCAAAAACTGAATCCCAAGACCAAACATTGTTATGCTGTGCTTCTTAAAATAAGTTATGAGATGGGAAAGAGGGCACCCCAAATATATATATATAATTATATATAATATAATATATATAATATACATAATTATGTATAATATAATATATATAATATACATAATCTCCTTTTACATCCTGCATCCTTATATTATATATAATATATTATATATAATAATATATAATAATATATAATTATTATATATATAATATATACAATTATATATGATTATATAATTGTATATAATTATATATATAATATAATAATATATAATATCTATTATTATATTATATATAATATAATATTATAATATATATCATATAT", "ctgggactacaccagctgccacatgcctggctaatttttcgtatttttagtagagacagggtttcactgaattggccaggatggtcttgatctcctcaccttgtgatcttctgccttggcctcccaaagtgctgggattacaggcctgagccaagtacatatttttaaaggaagaaaaatttcaaaggtactctgtttggtacaataatcaaatatataaattgaggaataaaacataaccacggaacatatttataactgcttatggaaaatacagaggataattttttaaataacatattttgaaaaccatttatttaggaatttgaaaagatcgcatttgacaggccagtatgaacacaacttgaacgcagcaagacaggttccccataaaaaaaaatcaaacacgtaagggaaaatgaaaccacaaagtttcaatctgctctgacctttgaaaactcagcacagacagtggcacttaggaccaacagcaggagatccctaatcccatcaccgtggcgatagggcataaagattccagggtgaagtcacaatccacattgggaggtccaactgctgccaggcagacaggtgtgcctttacatgtacaggaaggtcattgaaggctcagtgttttgtttcaaaaactgaatcccaagaccaaacattgttatgctgtgcttcttaaaataagttatgaggaatgggaaagagggcaccccaaatatatatatataattatatataatataatatatataatatacataattattataatataaatatatataatatacataatctccttttacatcctgcatccttatattatatataatatattatatataataatatataataatatataattattatatatataatatatacaattatataagattatataattgtatataattatatatataatataataatatataatagctattattatattatatataatataatatttaatatatatcat", 0, 100, 18);
  }     //                          CTGGGACTACACCAGCTGCCACATGCCTGGCTAATTTTTCGTATTTTTAGTAGAGACAGGGTTTCACTGAATTGGCCAGGATGGTCTTGATCTCCTCACCTTGTGATCTTCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCCTGAGCCAAGTACATATTTTTAAAGGAAGAAAAATTTCAAAGGTACTCTGTTTGGTACAATAATCAAATATATAAATTGAGGAATAAAACATAACCACGGAACATATTTATAACTGCTTATGGAAAATACAGAGGATAATTTTTTAAATAACATATTTTGAAAACCATTTATTTAGGAATTTGAAAAGATCGCATTTGACAGGCCAGTATGAACACAACTTGAACGCAGCAAGACAGGTTCCCCATAAAAAAAAATCAAACACGTAAGGGAAAATGAAACCACAAAGTTTCAATCTGCTCTGACCTTTGAAAACTCAGCACAGACAGTGGCACTTAGGACCAACAGCAGGAGATCCCTAATCCCATCACCGTGGCGATAGGGCATAAAGATTCCAGGGTGAAGTCACAATCCACATTGGGAGGTCCAACTGCTGCCAGGCAGACAGGTGTGCCTTTACATGTACAGGAAGGTCATTGAAGGCTCAGTGTTTTGTTTCAAAAACTGAATCCCAAGACCAAACATTGTTATGCTGTGCTTCTTAAAATAAGTTATGAGGAATGGGAAAGAGGGCACCCCAAATATATATATATAATTATATATAATATAATATATATAATATACATAATTATTATAATATAAATATATATAATATACATAATCTCCTTTTACATCCTGCATCCTTATATTATATATAATATATTATATATAATAATATATAATAATATATAATTATTATATATATAATATATACAATTATATAAGATTATATAATTGTATATAATTATATATATAATATAATAATATATAATAGCTATTATTATATTATATATAATATAATATTTAATATATATCAT


  public void checkEarlyTerminationDoesntLie(String templateStr, String readStr, int startPos, int maxScore, int maxShift) {
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString(templateStr);
    final boolean rc = false;
    final int length = read.length;
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    final EditDistance e = EditDistanceFactory.createEditDistance(params, null, length, length);

    int[] actions = e.calculateEditDistance(read, length, template, startPos, rc, maxScore, maxShift, true);
    final int initialAS = actions[ActionsHelper.ALIGNMENT_SCORE_INDEX];
    // if (initialAS < maxScore) {
    //   System.err.println("initial AS " + initialAS + " so alignment with score < " + maxScore + " possible at pos " + actions[ActionsHelper.TEMPLATE_START_INDEX]);
    // } else {
    //   System.err.println("initial AS " + initialAS + " so no low scoring alignment possible");
    // }

    // final LowerBoundEditDistance le = new LowerBoundEditDistance(5);
    // final int[] actions2 = le.calculateEditDistance(read, length, template, start, 100);
    // final int newAS2 = actions2[ActionsHelper.ALIGNMENT_SCORE_INDEX];
    // System.err.println("LB AS " + newAS2 + " at pos " + actions2[ActionsHelper.TEMPLATE_START_INDEX]);


    // Try finding alignment with a looser lower bound
    actions = e.calculateEditDistance(read, length, template, startPos, rc, 1000000, maxShift, true);
    final int newAS = actions[ActionsHelper.ALIGNMENT_SCORE_INDEX];
    //System.err.println(ActionsHelper.toString(actions));
    //System.err.println("Actual alignment AS " + newAS + " at pos " + actions[ActionsHelper.TEMPLATE_START_INDEX]);
    if (initialAS >= maxScore) {
      assertTrue("Terminated early with maxscore " + maxScore + " but it actually could find an alignment with score " + newAS + "!", newAS >= maxScore);
    } else {
      assertTrue("Gave different alignment scores depending on max score score", initialAS == newAS);
    }
  }

  private static final String CHR22_TEMPLATE = "TTCTCAGTGACTTATGGTAAATGACTTAGCAGCTTTTTTTTTTTTTTTTTTTTTTTCTGTGAAGAGTCTTCCTCTGTTGCCCAGGCTGGAGTGCAGTTGGCACAATTTCGGCTCATTGTGACCTCTGCCTCCTGGGCTTAAGTGATCCTCCCACCTCAGCTTCCCAGTAGCTGGGACTACAGGTGCGCACCACCACACGCGGCTTTTTTTTTTCTTCTTTTTTTTGTAGAGATCAAGTTTCACTATATTACCCAGGCTGATCTCCAACTGGGCTCCTCTGCCTCCCAGAGTGCTGGGACAGGGTCTCACTCTGTTGCCCAGGCTGGAGTGCAGTGGCATCGTCCCGGCTCACTGCAGCCTAGACCTCCCAGGCTCAAGTGATCCTCCTGCCTCAGCCTCCTGAGTAGCTGGAACCACATGCGTGCACCACCACACGCAGCGAACTATTTTGTA";
  private static final String READ_FOO = "CCAACTGCACTCCAGCCTGGGAGCATAAGACCCTGTCAAAAAGGAAAACAAAAAAAAAAATAAAAAAGCCATGTCAAAGCCAAGAAGCAGGATATTCTCT";

  public void testMaxShift() {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 100, 100);
    final byte[] read = DnaUtils.encodeString(READ_FOO);
    final byte[] template = DnaUtils.encodeString(CHR22_TEMPLATE);
    final int[] actions = ed.calculateEditDistance(read, 100, template, 0, true, 10, 8, true);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }


  public void testMaxShift2() {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 100, 100);
    final byte[] read = DnaUtils.encodeString("GCCCAGTTAATTTTTTCTATTTTTAGTTGAGACTGGGTTTCACCACGTTGGCCAGGCTCCTGACCTCAGGTGAGTCACCTGCCTTGGCCTCCCAAAGTGC");
    final byte[] template = DnaUtils.encodeString("GCCCAGTTAATTTTTTCTATTTTTAGTTGAGACTGGGTTTCACCATGTTGGCCAGGCTGGGGAATTCCTGACCTCAGGTGAGTCACCTGCCTTGGCCTCCCAAAGTGCTGGGGTTACAAGCATGAGTCACCGCACCTGGCCACTTCCCACAG");
    int[] actions = ed.calculateEditDistance(read, 100, template, 0, false, 20, 8, true);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(10, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

    actions = ed.calculateEditDistance(read, 100, template, 0, false, 20, 7, true);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  public void testCreateProtein() throws InvalidParamsException, IOException {
    assertNotNull(EditDistanceFactory.createProteinEditDistance(new ProteinScoringMatrix()));
  }

  public void testMNPAlignment() {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 100, 100);

//    596     83      gi|89161203|ref|NC_000022.9|NC_000022   139     255     95=8I18=8D79=   =       2284146 -206    ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:18 NM:i:16 MQ:i:255        XA:i:37 IH:i:1
//    224561  99      gi|89161203|ref|NC_000022.9|NC_000022   147     255     91=5X99=1X4=    =       2284178 218     GAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCACTGCTAGCT        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:6  NM:i:6  MQ:i:255        XA:i:24 IH:i:1

    final String readStr =                                "ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC";
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString("TTGGTGTCATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCACTGCAAGCTCCGCTTCCC");

    int[] actions = ed.calculateEditDistance(read, read.length, template, 8, false, 30, 8, false);
    assertEquals(8, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(5, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

    final byte[] readRc = DnaUtils.encodeString(DnaUtils.reverseComplement(readStr));
    actions = ed.calculateEditDistance(readRc, read.length, template, 8, true, 30, 8, false);
    assertEquals(8, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(5, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    ed.logStats();

  }


  public void testOffTemplateAlignment() {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                     .create();
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 101, 101);

//    596     83      gi|89161203|ref|NC_000022.9|NC_000022   139     255     95=8I18=8D79=   =       2284146 -206    ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:18 NM:i:16 MQ:i:255        XA:i:37 IH:i:1
//    224561  99      gi|89161203|ref|NC_000022.9|NC_000022   147     255     91=5X99=1X4=    =       2284178 218     GAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCACTGCTAGCT        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:6  NM:i:6  MQ:i:255        XA:i:24 IH:i:1

    final String readStr = "AGGTTTTTGGTTTTTTTGGTTTTTATTTTTTTGAGATGGAGTCTTGCACTGTCGCCTGGGCTGGAGTGCAATGGTGCGATCTTGGCTCACTGCAACCTCTG";
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString("NNNNNNNNNNNNNNNNNNNNCCGGGAGGTGGAGGTTGCAGTGAGCCGAGATCGCACCATTGCACTCCAGCCCAGGCGTCAGGGCAAGACTCCATCTCAAATAAAGAAGAAGAAAATGAAACNNNNNNNNNNNNNNNNNNNN");

    final int[] actions = ed.calculateEditDistance(read, read.length, template, 28, true, 2, 8, false);
    assertEquals(28, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

  }

//  public void testOffTemplateAlignment2() {
//
//    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(19).gapExtendPenalty(1).substitutionPenalty(9).nsAsMismatches(true)
//                                     .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
//                                     .create();
//    final String readStr = "ATTAAATTGGTACCTTGTGAATTCTGATAGCATGAAATAACATGCTATATTGAAAAGTTGCATGTTAACAGACTGAGACTGC";
//    final EditDistance ed = new RcEditDistance(new SingleIndelEditDistance(params, readStr.length())); // EditDistanceFactory.createEditDistance(params, null, readStr.length(), readStr.length());
//
//    final byte[] read = DnaUtils.encodeString(readStr);
//    final byte[] template = DnaUtils.encodeString("CTCAGTCTGTTAACATGCAACTTTTCAATATAGCATGTTATTTCATGCTATCAGAATTCACAAGGTACCAATTTAATTACTACAGAGTACTTATAGAATCATTTAAAATATAATAAAATTGTATGATAGAGATTATATGCAATAAAACATTAACAAAATGCTAAAATACGAGACATATTCCGATTAAAGTATTTATAAAA");
//
//      final int[] actions = ed.calculateEditDistance(read, read.length, template, 0, true, 50, 14, false);
//
//  }

  public void testAlignerChainAutoIllumina() {
    // Should end up using the single indel aligner, which will have heaps of mismatches near the end.
    final SAMReadGroupRecord readgroup = new SAMReadGroupRecord("foo");
    readgroup.setPlatform("ILLUMINA");
    final AlignerMode alignerMode = AlignerMode.AUTO;
    final NgsParams params = getNgsParams(readgroup, alignerMode);

    final int[] actions = alignerChainTestRead(params);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(110, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  private NgsParams getNgsParams(SAMReadGroupRecord readgroup, AlignerMode alignerMode) {
    return new NgsParamsBuilder().gapOpenPenalty(5).gapExtendPenalty(1).substitutionPenalty(15).alignerMode(alignerMode)
          .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).readGroup(readgroup).create())
          .create();
  }

  public void testAlignerChainGeneralIlluminaOverride() {
    // Should use the general purpose aligner, two inserts
    final SAMReadGroupRecord readgroup = new SAMReadGroupRecord("foo");
    readgroup.setPlatform("ILLUMINA");
    final NgsParams params = getNgsParams(readgroup, AlignerMode.GENERAL);

    final int[] actions = alignerChainTestRead(params);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(18, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  public void testAlignerChainAutoNotIllumina() {
    // Auto aligner selection shouldn't pick the illumina one
    final SAMReadGroupRecord readgroup = new SAMReadGroupRecord("foo");
    readgroup.setPlatform("NotIllumina");
    final NgsParams params = getNgsParams(readgroup, AlignerMode.AUTO);

    final int[] actions = alignerChainTestRead(params);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals(18, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  private int[] alignerChainTestRead(NgsParams params) {
    // Align a read with two inserts with the aligner determined by the params
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 10, 101);
    final String readStr = "GAGGTTGTTTTCAGTGAGCCGAGATCGCACCATTGCACTCCAGCCCAGGCGTCAGGGCAAGACTCCATCTCATTTTAATAAAGAAGAAG";
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString("CCGGGAGGTGGAGGTTGCAGTGAGCCGAGATCGCACCATTGCACTCCAGCCCAGGCGTCAGGGCAAGACTCCATCTCAAATAAAGAAGAAGAAAATGAAAC");

    return ed.calculateEditDistance(read, read.length, template, 10, false, 200, 8, false);
  }

  public void testWoogieStartPosition() {
    final NgsParams params = new NgsParamsBuilder().singleIndelPenalties("alignertable").alignerMode(AlignerMode.TABLE)
        .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
        .create();
    final EditDistance ed = EditDistanceFactory.createEditDistance(params, null, 101, 101);

    final String readStr = "GCATGGACTGGTTTGTAGATGATTCGACTTTAATTGGCCCGAGGTCTATACGCGATGGCTTTTGTACTCCTCTGGCGATGAAGCCGTCTTGGTACCGTATT";
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString("GGGAAGTGCGTACAACGTATGGGTCTGAATAGGTGTTTGTCTGGGCTATAGCTTATAAGCAGCTTGGATTTATAGACATTTAGAGCTTTTCATTAGTCGGGCATGGACTGGTTTGTAGATGATTCGACTTTAATTGGCCCGAGGTCTATACGCGGTGGCTTTTGTACTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

    // At one point we had problems with results that returned Integer.MAX_VALUE score but reset the template start position to 0.
    final int[] actions = ed.calculateEditDistance(read, read.length, template, 100, false, 90, 50, false);
    assertNotNull(actions);
    assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(100, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }
}
