/*
 * Copyright (c) 2016. Real Time Genomics Limited.
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

import static com.rtg.mode.DnaUtils.reverseComplement;

import java.io.IOException;

import org.junit.Assert;

import com.rtg.alignment.SingleIndelSeededEditDistance;
import com.rtg.alignment.UnidirectionalAdaptor;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.ngs.MapParamsHelper;
import com.rtg.ngs.NgsParams;
import com.rtg.util.cli.CFlags;

/**
 */
public class PairAlignerTest extends AbstractNanoTest {

  public void testAlignmentOverlappingNoReadThrough() {
    final PairAligner pairAligner = getPairAligner(0);
    final String leftString =  "ACGATGACGTAGATGTGTACCGTTATATGC";
    final String rightString = "----TGACGTAGATGTGTACCGTTATATGCCCT";
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", leftString), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(rightString)), right.toFasta());
  }
  private String fasta(String name, String read) {
    return ">" + name + "\n" + read.replaceAll("-", "") + "\n";
  }

  public void testAlignmentTrimSingleBaseFromR1() {
    final PairAligner pairAligner = getPairAligner(0);
    final String leftString =  "ACGATGACGTAGATGTGTACCGTTATATGCCCTG";
    final String rightString = "----TGACGTAGATGTGTACCGTTATATGCCCT";
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", leftString.substring(0, leftString.length() - 1)), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(rightString)), right.toFasta());
  }

  public void testAlignmentTrimSingleBaseFromR2() {
    final PairAligner pairAligner = getPairAligner(0);
    final String leftString =  "-ACGATGACGTAGATGTGTACCGTTATATGCC";
    final String rightString = "GACGATGACGTAGATGTGTACCGTTATATGCCCT";
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", leftString), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(rightString.substring(1))), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughBothEnds() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "GCCG";
    final String middle = "ACGATGACGTAGATGTGTACCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = middle + suffix;
    final String rightString = prefix + middle;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middle)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughBothEndsLongerR2() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "CGCCG";
    final String middle = "ACGATGACGTAGATGTGTACCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = middle + suffix;
    final String rightString = prefix + middle;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middle)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughBothEndsShorterR2() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "CGCCG";
    final String middle = "ACGATGACGTAGATGTGTACCGTTATA";
    final String suffix = "TGCAGAT";
    final String leftString = middle + suffix;
    final String rightString = prefix + middle;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middle)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughR1LongerR1() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "ACCGCCG";
    final String middle = "ACGATGACGTAGATGTGTACCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = prefix + middle + suffix;
    final String rightString = middle;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    // Should trim just the suffix off the left read
    Assert.assertEquals(fasta("left", prefix + middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middle)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughR2LongerR2() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "ACCGCCG";
    final String middle = "ACGATGACGTAGATGTGTACCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = middle;
    final String rightString = prefix + middle + suffix;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    // Should trim just the prefix off the right read
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middle + suffix)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughBothEndsMismatches() {
    final PairAligner pairAligner = getPairAligner(0, PairAligner.MismatchType.PICK_BEST);
    final String prefix = "GCCG";
    final String middle =          "ACGATGACGTAGA-T-GTGT-A-CCGTTATATGC";
    final String middleMismatch =  "ACGATGACGTAGA-C-GTGT-T-CCGTTATATGC";
    final String middleMismatch2 = "ACGATGACGTAGA-T-GTGT-A-CCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = middle + suffix;
    final String rightString = prefix + middleMismatch;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString, 40);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString), 30);
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middleMismatch2)), right.toFasta());
  }

  public void testAlignmentOverlappingReadThroughBothEndsInsertsInR2() {
    final PairAligner pairAligner = getPairAligner(0);
    final String prefix = "GCCG";
    final String middle =       "ACGATGACGTAGATGT--GTACCGTTATATGC";
    final String middleInsert = "ACGATGACGTAGATGTGGGTACCGTTATATGC";
    final String suffix = "AGAT";
    final String leftString = middle + suffix;
    final String rightString = prefix + middleInsert;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middleInsert)), right.toFasta());
  }

  private PairAligner getPairAligner(int probeLength) {
    return getPairAligner(probeLength, PairAligner.MismatchType.NONE);
  }

  private PairAligner getPairAligner(int probeLength, PairAligner.MismatchType t) {
    final int maxReadLength = 300;
    final int seedLength = 5;
    final CFlags flags = new CFlags();
    PairedEndTrimCli.initFlags(flags);
    final NgsParams ngsParams = MapParamsHelper.populateAlignerPenaltiesParams(NgsParams.builder(), flags)
      .singleIndelPenalties(null)
      .create();
    return new PairAligner(
      new UnidirectionalAdaptor(new SingleIndelSeededEditDistance(ngsParams, false, seedLength, 80, 80, maxReadLength)),
      25, 90, probeLength, 0, 0, false, t, false);
  }


  public void testRealWorld() {
    // Observed in MI000843
    final PairAligner pairAligner = getPairAligner(0);
    final String leftString = "AAATACTCAGGACATGTCAGCCTCTCAGGTTGATGTAGCTGTGAAAATTAATAAGAAAGTTGTGCCCCTGGACTTTTCTATGAGTTCTTTAGCTAAACGAATAAAGCTGTT";
    final String rightString = "CCAGAAGGGACTCCAGCTCAGACGTGTGCTCTTCCGATCTAAATACTCAGGACATGTCAGCCTCTCAGGTTGATGTAGCTGTGAAAATTAATAAGAAAGTTGTGCCCCTGGACTTTTCTATGAGTTCTTTAGCTAAACGAATAAAGCAGTT";
    final String expectedLeft = leftString;
    // Expect to trim first 40 bases of right read
    final String expectedRight = rightString.substring(40);
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", expectedLeft), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(expectedRight)), right.toFasta());
  }

  public void testRealWorld2() {
    // Observed in MI000843
    final PairAligner pairAligner = getPairAligner(0);
    final String leftString = "ATGATAGATTTATCGCTTCTGTGACAGACAGTGAAAACACAAATCAAAGAGAAGCTGCAAGTCATGGTAAGTCCTCTGTTTAGTTGAACTACAGGTTTTTTTGTTGTTGTTGT";
    final String rightString =
      "ACAAAAACAACAAAAAAACCTGTAGTTAAACAAAAAAGAGGACGTAACATGAATTGGAGCGTGTAATTTAATTTTGTTTTAACAGTCTGTAACAGAAACGATAAAAAAATAAAAGAAAGGAAGAGCACAAGTATGAAAAACAATAACACA";
    final String expectedLeft = leftString;
    final String expectedRight = rightString;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", rightString);
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", expectedLeft), left.toFasta());
    Assert.assertEquals(fasta("right", expectedRight), right.toFasta());
  }

  public void testRealProbeOverlap() {
    // Observed in MI000843
    final PairAligner pairAligner = getPairAligner(37);
    final String leftString = "GTCTGTGGGATCTGGTACAGTTTTGTAAATTCATCTGGATAATTTTCTACCCAGTTCCAAAATGCCTATAAATAAAGAACCACATGTCTTTACAGTGATTTTTTGTAATACAGTAAATAAAATGTCACAA";
    final String rightString = "AAGTTTGTGACATTTTATTTACTGTATTACAAAAAATCACTGTAAAGACATGTGGTTCTTTATTTATAGGCATTTTGGAACTGGGTAGAAAATTATCCAGATGAATTTACAAAACTGTACCAGATCCCAC";
    final String expectedLeft = leftString;
    final String expectedRight = rightString.substring(0, rightString.length() - 33);
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", rightString);
    pairAligner.processReads(new FastqPair(left, right));
    Assert.assertEquals(fasta("left", expectedLeft), left.toFasta());
    Assert.assertEquals(fasta("right", expectedRight), right.toFasta());
  }

  public void testDeletionZeroPhred() throws IOException {
    final PairAligner pairAligner = getPairAligner(0, PairAligner.MismatchType.ZERO_PHRED);
    //                          0123456789012345678901234567890123456
    final String leftString =  "ATGATAGATAAATCAAAGAAGCTGCTATAGGCGGTGGCTACC";
    //                          012345678--------90123456
    final String rightString = "ATGATAGATAA----AAGAAGCTGCTATAGGCGGTGGCTACC";
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString, leftString.length());
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString), leftString.length());

    // Make sure this is treated as an deletion
    Assert.assertTrue(left.length() > right.length());
    pairAligner.processReads(new FastqPair(left, right));
    mNano.check("pair-aligner-zero-phred-del.fq", left.toFastq() + right.toFastq());
    Assert.assertEquals(fasta("left", leftString), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(rightString)), right.toFasta());
  }

  public void testInsertionZeroPhred() throws IOException {
    final PairAligner pairAligner = getPairAligner(0, PairAligner.MismatchType.ZERO_PHRED);
    final String prefix = "";
    final String middle =       "ACGATGACGTAGATGT--GTACCGTTATATGC";
    final String middleInsert = "ACGATGACGTAGATGTGGGTACCGTTATATGC";
    final String suffix = "";
    final String leftString = middle + suffix;
    final String rightString = prefix + middleInsert;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString, leftString.length());
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString), rightString.length());
    pairAligner.processReads(new FastqPair(left, right));
    mNano.check("pair-aligner-zero-phred-ins.fq", left.toFastq() + right.toFastq());
    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middleInsert)), right.toFasta());
  }


  public void testInsertion() {
    final PairAligner pairAligner = getPairAligner(0);
    final String middle = "TGAGATTACGAT---GATAGATAAATCAAA";
    final String suffix = "GAAGCTGC";
    final String leftString =  middle + suffix;
    final String prefix = "ACC";
    final String middleInserted = "TGAGATTACGATGGGGATAGATAAATCAAA";
    final String rightString = prefix + middleInserted;
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", reverseComplement(rightString));

    // Make sure this is treated as an insertion
    Assert.assertTrue(left.length() > right.length());

    pairAligner.processReads(new FastqPair(left, right));

    Assert.assertEquals(fasta("left", middle), left.toFasta());
    Assert.assertEquals(fasta("right", reverseComplement(middleInserted)), right.toFasta());
  }

  public void testFullyTrimmedR1() {
    final PairAligner pairAligner = getPairAligner(30);
    final String leftString = "GCCTGTAACCTAGAAATGGGACAGA";
    final String rightString = "GTCCCATTTCTAGGTTACAGGCAGA";
    final FastqSequence left = FastqSequenceTest.getFastq("left", leftString);
    final FastqSequence right = FastqSequenceTest.getFastq("right", rightString);

    // Make sure this is treated as an insertion
    pairAligner.processReads(new FastqPair(left, right));

    Assert.assertEquals(fasta("left", leftString.substring(0, leftString.length() - 3)), left.toFasta());
    Assert.assertEquals(fasta("right", ""), right.toFasta());
  }
}
