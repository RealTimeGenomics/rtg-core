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

import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.protein.GotohProteinEditDistanceTest;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;


/**
 */
public class ActionsValidatorTest extends TestCase {

  public void testNull() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    assertFalse(av.isValid(null, null, 0, null, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(null, null, 0, null, 0);
    assertEquals("ValidatingEditDistance action problem: actions array is null" + StringUtils.LS, error);
  }

  public void testTooShort() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    assertFalse(av.isValid(new int[2], null, 0, null, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(new int[2], null, 0, null, 0);
    assertEquals("ValidatingEditDistance action problem: actions array is too short: 2" + StringUtils.LS, error);
  }

  public void testNumberOfActions() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final int[] actions = new int[3];
    actions[ActionsHelper.ACTIONS_LENGTH_INDEX] = 4;
    assertFalse(av.isValid(actions, null, 0, null, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(actions, null, 0, null, 0);
    assertEquals("ValidatingEditDistance action problem: int[] is shorter than number of actions requires" + StringUtils.LS, error);
  }

  public void testMaxIntAction() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final int[] actions = new int[3];
    actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
    assertTrue(av.isValid(actions, null, 0, null, 0));
  }

  public void testSameFail() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aata");
    final int[] actions = ActionsHelper.build("====", 0, 0);
    assertFalse(av.isValid(actions, read, read.length, tmpl, Integer.MAX_VALUE));
    String error = av.getErrorDetails(actions, read, read.length, tmpl, 0);
    TestUtils.containsAll(error, "ValidatingEditDistance action problem: action 3: read[2] (1) != template[2] (4) score=0",
        "0|       10|       20|       30|       40|       50|       60|       70|       80|       90|",
        "tmpl:    AATA", "read:    AAAA", "actions: ==== score=0");

    error = av.getErrorDetails(actions, read, read.length, tmpl, 1);
    TestUtils.containsAll(error, "tmpl[0..]: A", "tmpl:    ATAN", "read:    AAAA", "actions: ==== score=0");
  }

  public void testEndInsert() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aata");
    final int[] actions = ActionsHelper.build("====", 0, 0);
    assertFalse(av.isValid(actions, read, 0, tmpl, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(actions, read, 0, tmpl, 0);
    assertTrue(error.contains("ValidatingEditDistance action problem: non-insert action at end of the read"));
  }

  public void testSubstitutionFail() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aaaa");
    final int[] actions = ActionsHelper.build("==X=", 0, 0);
    assertFalse(av.isValid(actions, read, read.length, tmpl, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(actions, read, read.length, tmpl, 0);
    assertTrue(error.contains("ValidatingEditDistance action problem: action 3: read[2] (1) == template[2] (1) score=1"));
  }

  public void testFewActionsFail() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aaaa");
    final int[] actions = ActionsHelper.build("===", 0, 0);
    assertFalse(av.isValid(actions, read, read.length, tmpl, Integer.MAX_VALUE));
    final String error = av.getErrorDetails(actions, read, read.length, tmpl, 0);
    assertTrue(error.contains("ValidatingEditDistance action problem: actions cover only 3 residues, but the read has 4"));
  }

  public void testClaimedScoreMaxIntFail() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aaaa");
    final int[] actions = ActionsHelper.build("====", 0, Integer.MAX_VALUE);
    assertFalse(av.isValid(actions, read, read.length, tmpl, 63));
    final String error = av.getErrorDetails(actions, read, read.length, tmpl, 0);
    assertTrue(error.contains("ValidatingEditDistance action problem: actual score 0 < max (63) but score was MAX_VALUE"));
  }

  public void testClaimedScoreFail() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aaaa");
    final int[] actions = ActionsHelper.build("====", 0, 5);
    assertFalse(av.isValid(actions, read, read.length, tmpl, 63));
    final String error = av.getErrorDetails(actions, read, read.length, tmpl, 0);
    assertTrue(error.contains("ValidatingEditDistance action problem: actual score 0 != claimed score 5"));
  }

  public void testClaimedScorePass() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("aaaa");
    final int[] actions = ActionsHelper.build("====", 0, 0);
    assertTrue(av.isValid(actions, read, read.length, tmpl, 63));
  }

  public void testClaimedScoreNPass() {
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 2);
    final byte[] read = DnaUtils.encodeString("aaaa");
    final byte[] tmpl = DnaUtils.encodeString("agNa");
    final int[] actions = ActionsHelper.build("=XX=", 0, 3);
    assertTrue(av.isValid(actions, read, read.length, tmpl, 63));
  }

  public void testEqual() {
    checkIsValidDna(
        "acgt",
        "acgt");
  }

  /**
   * Test a deletion from the read.
   */
  public void testDelete() {
    checkIsValidDna(
        "acgacgtttcgcgcgc",
        "cgacgcgcgcg");
    }

  /**
   * Two deletes from the read
   */
  public void testDelete2() {
    checkIsValidDna(
        "tactcgaacccttcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaatatgggtactgcat",
        "tactcga    ttcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaatat   tactgcat".replaceAll(" ", ""));
  }

  /**
   * Test a single insertion into the read.
   */
  public void testInsert() {
    checkIsValidDna(
        "gatacgaactcgtacgcact",
        "gatacgaacccctcgtacgcactcg");
  }

  /**
   * This is an example where CachedMatrixEditDistance gives an alignment
   * with score 5, and misses the <code>8=3D11=</code> alignment with score 4.
   */
  public void testBadAlignment() {
    checkIsValidDna(
        "gatacgaa" +    "ctcgtcgcact",
        "gatacgaa" + "accctcgtcgcactcg");
  }

  /**
   * Test two insertions.
   */
  public void testInsert2() {
    checkIsValidDna(
           "actcactaagggggtttctatagtttttcactcgg",
        "gggactcactaa" + "tttctatag" + "cactcggggg");
  }

  public void testProteinEquals() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHILKMFPSTWYV",
        "ARNDCQEGHILKMFPSTWYV"
        );
  }

  public void testProteinSubs() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHILKMFPSTWYV",
        "ARNACQEGHLIKMGPSTWYV"
        );
  }

  public void testProteinDelete() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQEGHIIILLLLKMFPSTWYV",
        "ARNACQEGHI" + "LKMGPSTWYV"
        );
  }

  public void testProteinInserts() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "ARNDCQE    GHILKMF    PSTWYV".replaceAll(" ", ""),
         "RNACQEAVAVGHILKMFFFFFPSTWY"
        );
  }

  /** An example where ProteinEditDistance returns a non-optimal alignment. */
  public void testProteinBadAlignment() throws InvalidParamsException, IOException {
    checkIsValidProtein(
        "CEQGIHIILLLLKMFPSTWYV",
        "CQEGHI" + "LKMGPSTWYV"
        );
  }

  public void testRcSimple() {
    final String read = "tactcgaaccccttctatggggtactgcat";
    checkIsValid(
        read,
        DnaUtils.reverseComplement(read),
        true,
      1);
  }

  public void testRcWithDelete() {
    checkIsValid(
        "acgacgtttcgcgcgc",
        DnaUtils.reverseComplement("cgacgcgcgcg"),
        true,
      1
        );
  }

  /*
   This test does not work yet - the calculation of RC start position is 2, but should be 0?
   public void testRcWithInsert() {
    checkIsValid(
        "gatacgaactcgtacgcact",
        Utils.reverseComplement("gatacgaacccctcgtacgcactcg"),
        true,
        Integer.valueOf(1)
        );
  }
  */

  public void checkIsValidDna(final String readString, final String templateString) {
    checkIsValid(readString, templateString, false, 1);
  }

  public void checkIsValidProtein(final String readString, final String templateString) throws InvalidParamsException, IOException {
    checkIsValid(readString, templateString, false, new ProteinScoringMatrix());
  }

  /**
   * Test the ActionsValidator by mutating entries in the actions array.
   *
   * @param readString the read
   * @param templateString the template
   * @param scoring a ProteinScoringMatrix, or an Integer (for the DNA gap open penalty).
   */
  public void checkIsValid(final String readString, final String templateString, final boolean rc, final Object scoring) {
    final byte[] read;
    final byte[] temp;
    final int maxScore = Integer.MAX_VALUE;
    final EditDistance editDist;
    final ActionsValidator validator;
    final MaxShiftFactor msf;
    if (scoring instanceof Integer) {
      final int gapOpen = (Integer) scoring;
      final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpen).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(1).create();
      editDist = new RcEditDistance(new GotohEditDistance(params));
      validator = new ActionsValidator(gapOpen, 1, 1, 1);
      // encode as DNA
      read = readString.getBytes();
      temp = templateString.getBytes();
      DnaUtils.encodeArray(read);
      DnaUtils.encodeArray(temp);
      msf = params.alignerBandWidthFactor();
    } else {
      final ProteinScoringMatrix matrix = (ProteinScoringMatrix) scoring;
      editDist = EditDistanceFactory.createProteinEditDistance(matrix);
      validator = new ActionsValidator(matrix);
      // encode as Protein
      read = GotohProteinEditDistanceTest.encodeProteins(readString);
      temp = GotohProteinEditDistanceTest.encodeProteins(templateString);
      msf = new MaxShiftFactor(0.7);
    }
    int[] a0 = editDist.calculateEditDistance(read, read.length, temp, 0, rc, maxScore, msf.calculateMaxShift(read.length), true);
    //System.out.println("score=" + a0[ActionsHelper.ALIGNMENT_SCORE_INDEX] + " alignment=" + ActionsHelper.toString(a0));
    a0 = ActionsHelper.copy(a0); // make it minimal.
    //Note: Only print alignment results if DNA.
    //AlignmentResult ares = new AlignmentResult(read, a0, 0, temp);
    //System.out.println(ares.toString());
    assertTrue(validator.isValid(a0, read, read.length, temp, rc, maxScore));

    // Now try perturbing each integer in the actions array,
    // and check that every change gives a false result from isValid.
    final int len = ActionsHelper.actionsCount(a0);
    for (int i = 0; i <= ActionsHelper.ACTIONS_START_INDEX + (len - 1) / ActionsHelper.ACTIONS_PER_INT; i++) {
      //System.err.println("\nPerturbing entry a0[" + i + "] = " + a0[i]);
      a0[i]--;
      assertFalse(validator.mErrorMsg, validator.isValid(a0, read, read.length, temp, rc, maxScore));
      a0[i]++;
      assertTrue(validator.mErrorMsg, validator.isValid(a0, read, read.length, temp, rc, maxScore));
      a0[i]++;
      assertFalse(validator.isValid(a0, read, read.length, temp, rc, maxScore));
      a0[i]--;
      assertTrue(validator.isValid(a0, read, read.length, temp, rc, maxScore));
    }
  }
}
